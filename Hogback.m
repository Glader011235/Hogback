%11/8/2017
%This code evolves a hogback through time. A resistant layer of
%rock, which weathers slowly, overlies a softer layer of rock that weathers
%quickly. Resistant rock produces "blocks" which land on the adjoining
%hillslope. Boundaries incise at a specified rate. User can set hogback
%layer thickness, block size, and dip, as well as relative weathering and
%incision rates. Trackable metrics included are time and space-averaged
%slope, block height, weathering rate, and erosion rate. Parameter space
%exploration allows use of different weathering rules and block movement
%rules. Parameters that users need to specify are surrounded by many 
%comment signs. Model parameters are set to match those presented in JGR
%paper (2017?). This model differs from that of Glade et al., 2017 in that
%it allows the boundary condition to move with the feature. 

%R.C. Glade and R.S. Anderson


clear all %clear all current variables
close all %close all open figures

%% Initialize

%User input needed:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hogback characteristics
thickness = 1; %thickness of resistant layer [m]
blockheight = 1; %height of blocks [m]
blockwidth = 1; %width of blocks [m]
slopedegrees = 30; %dip of resistant layer [degrees]. Beware extremely low or high angles (must change initial topography setup to work)
drop = 3; %elevation drop (between resistant layer and slope) required for blocks to be released
blockmove = 2; %elevation difference between block and next cell required for block to move. E.g., if blockmove = 2, then elevation difference must be 2 times the current block height
cutoff = 0.2; %blocks no longer act as blocks after this block height
metric = 1; %toggle calculations of metrics in the model (increases run time). 0 = off; 1 = on.  
channel_on = 1; %=1 allows boundary to move with feature; =0 keeps boundary stationary

%Weathering parameters
exponential = 1; %=1 for exponential weathering, =0 for humped weathering rule
wstar = 0.1; %characteristic weathering depth
wdotnot = 1e-3; %maximum (bare bedrock) weathering rate
wdotnot_hard = 1e-5; %maximum (bare bedrock) weathering rate for resistant layer 
H0 = 0; %initial soil thickness

%for humped weathering curve:
%A = 5e-3; %peak of humped weathering curve
%A_l2 = A/100; %peak of humped weathering curve for resistant layer

%Transport parameters
k = 0.5; %this will be in units of [L/T] for Johnstone & Hilley style flux rule
hstar = 0.2; %characteristic depth for depth-dependent flux rule [m]
edot_channel = 3e-5; %channel boundary conditions at edges [m/yr]
%rockdensity = 2.3; %density of rock 
%soildensity = 2; %density of soil
%stdev = 5; %variable used for alternate spatial block distribution on slope

%Time conditions
dt = 2;
tmax = 1500000;
t = 0:dt:tmax;
imax = length(t);
nplots = 600;
tplot = tmax/nplots;
tchannel_stabilize = tmax;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%no user input needed (but there are some things you can change if desired)
%Block release rules 
hogbackslope = tand(slopedegrees); %convert slope degrees to percent
dx = blockwidth*cosd(slopedegrees); %set dx based on block height
elevhogdrop = blockheight*sind(slopedegrees); %calculate the vertical decrease in height on resistant layer after new row of blocks falls
xchange = dx;

%Domain setup
xmin = -250;
xmax = 150;
x = xmin:dx:xmax;
xedge = x(1:end-1)+(dx/2); %location of cell edges
zmax = 300; %peak elevation of initial topography
xlimfixed = [min(x) max(x)]; %for plotting purposes
ylimfixed = [123 320]; %for plotting purposes

%initial topography 
fractureplanedegrees = 90-slopedegrees; %line perpendicular to resistant dip slope to define "Fractures"
fractureplane = tand(fractureplanedegrees); %convert to percent slope
z1 = hogbackslope.*(x)+zmax; %line defining the top of the resistant layer
z2 = -fractureplane.*x+zmax; %line defining fracture plane
z2end = zmax-thickness*cosd(slopedegrees); %elevation of end of resistant layer depending on layer thickness
z3 = ones(1,length(x))*z2end; %horizontal line for initial topography of soft rock
zb = min(z1,z2); %combine to get initial topography
zb = max(zb,z3); %combine to get initial bedrock topography 
ztopo0 = zb; % saves initial topo


%Topography and lines for plotting 

%intercepts
b1 = zmax; %intercept for top of resistant layer
b2 = zmax-(thickness/(cosd(slopedegrees))); %intercept for bottom of resistant layer

% equations of lines bounding layers
z1l = (hogbackslope.*x) + b1; %line for top of resistant layer
z2l = z1-thickness/cosd(slopedegrees); %line for bottom of resistant layer
l2 = find(zb<=z1l & zb>=z2l); %location of the hard layer for indexing purposes



%Block rules 
blockend = length(x); %variable used later to prevent blocks going beyond domain space
nrocks = thickness/blockheight; %number of blocks released in one event
blockstart = l2(end)+1; %index of the beginning of the slope i.e. location of deposition of first block when released
wdotdist = []; %initialize array to track location of each block

%initial soil conditions
H=zeros(size(x));
H(:) = H0; %initial soil thickness

z = zb+H; %initial topography is sum of bedrock elevation and soil thickness

%other initializations
F(imax)=struct('cdata',[],'colormap',[]);
wdottracking=[]; %initialize array to keep track of block heights
peak=find(z==max(z(l2)));
elevblock=zeros(size(wdotdist));
elevnext=zeros(size(wdotdist));
Hsave=[]; %array to remember soil thickness underneath blocks


nframe = 0;
%slopelength=zeros(1,length(t)); %slope length tracking
%hogposition=zeros(1,length(t)); %peak position tracking


%plot initial topography
figure(1)
plot(x,zb)
axis equal

%initialize metric arrays
blockstart = 0; 
l=blockstart:5:blockstart+100; %used for calculating metrics
meanslope = zeros(1,length(l));
meanspacing = zeros(1,length(l));
meanH = zeros(1,length(l));
meanheight = zeros(1,length(l));
meanwdot = zeros(1,length(l));
meanflux = zeros(1,length(l));
erosion = zeros(1,length(l));
meanwdot2 = zeros(1,length(l));
analyticalslope = zeros(1,length(l));

meanslopesave=[];
analyticalslopesave=[];
meanfluxsave=[];
meanHsave=[];
meanheightsave=[];
meanwdotsave=[];
meandhdtsave=[];
erosionsave=[];
meanwdot2save=[];
blocklocationsave=[];
blockweathersave = [];
weathering_base_save = [];
peak_save = [];
meanerode_save = [];
meanslope2save = [];



channel = length(zb);





%% run

%Time loop

for i=1:imax
    

    %reevaluate location of hard layer
    l2=find(zb<=z1l & zb>=z2l); %this is the hard layer

    %set weathering rules:
    %if using humped weathering rule:
        if exponential==0

            wdot = wdotnot.*exp(-H./wstar)+((A*H./wstar).*exp(-H./wstar));
            wdot(l2) = 0; %don't allow resistant layer to weather- it only produces blocks. This conserves mass, but may be altered

        end

    %if using exponential weathering rule
        if exponential==1

            wdot = wdotnot.*exp(-H./wstar); 
            wdot(l2) = 0; %don't allow resistant layer to weather- it only produces blocks. This conserves mass, but may be altered

        end

        if channel_on == 1
            Hsave(find(wdotdist==min(channel)))=[]; %remove saved soil height from array at old block locations 
            wdottracking(find(wdotdist==min(channel)))=[]; %remove block weathering tracking from edge of domain
            wdotdist(find(wdotdist==min(channel)))=[]; %remove block location tracking from edge of domain
            wdottracking(find(wdotdist==min(channel)))=[]; %remove block weathering tracking from edge of domain
            wdotdist(find(wdotdist==min(channel)))=[]; %remove block location tracking from edge of domain

            Hsave(find(wdotdist==min(channel)-1))=[]; %remove saved soil height from array at old block locations 
            wdottracking(find(wdotdist==min(channel)-1))=[]; %remove block weathering tracking from edge of domain
            wdotdist(find(wdotdist==min(channel)-1))=[]; %remove block location tracking from edge of domain
            wdottracking(find(wdotdist==min(channel)-1))=[]; %remove block weathering tracking from edge of domain
            wdotdist(find(wdotdist==min(channel)-1))=[]; %remove block location tracking from edge of domain



        end
    

    %Moving blocks
    for j = 1:length(wdotdist)
        
        elevblock(j)=z(wdotdist(j)); %find elevation of top of block
        elevnext(j)=z(wdotdist(j)+1); %find elevation of cell downslope of block
        blockdrop=elevblock(j)-elevnext(j); %find elevation difference between top of block and next cell

            if blockdrop>=blockmove*(blockheight-wdottracking(j)) %if elevation difference is sufficient to move block   

                wdotdist(j)=wdotdist(j)+1; %moves the block one cell downslope in block location array
                H(wdotdist(j)-1)=Hsave(j); %restores previous soil height to cell where block used to be
                zb(wdotdist(j))=zb(wdotdist(j))+(blockheight-(wdottracking(j)))+H(wdotdist(j)); %Add block height to bedrock elevation of new block cell, including soil thickness
                zb(wdotdist(j)-1)=zb(wdotdist(j)-1)-(blockheight-(wdottracking(j)))-Hsave(j); %Subtract block height from bedrock elevation of previous block cell
                Hsave(j)=H(wdotdist(j)); %store soil thickness in new block cell
                H(wdotdist(j))=0; %set soil thickness to 0 on top of blocks
                z=H+zb; %update topography

            end
       
    end
   
    

    %Supply of blocks from hogback

    %keep track of elevation drop next to hogback
    elevdrop=z(l2(end))-z(l2(end)+1); %current elev diff between hogback and slope
    overshoot=elevdrop-drop; %how much further than one blockheight did slope decrease?

    if overshoot>=2*drop %check to make sure slope hasn't eroded down more than 1 block set. If so, need a smaller time step. 
    
        disp('DANGER')
    
    end

    peak=find(z==max(z(l2))); %find new peak elevation

    %block release
    if elevdrop>=drop & t(i)>0 %if elevation drop between resistant layer and slope is sufficient
        channel = [min(channel)-1, channel];
        zb(min(channel)) = zb(min(channel)+1);
        z(min(channel)) = z(min(channel)+1);
        

        
        zmax = z(peak)-elevhogdrop; %new peak after blocks fall based on block height
        zl2end=z(l2(end))-elevhogdrop; %new elev of end of resistant layer after blocks fall
        zb(peak-xchange/dx:l2(end)-xchange/dx)=linspace(zmax,zl2end,length(peak-xchange/dx:l2(end)-xchange/dx)); %place new zb after blocks fall
        zb(l2(end))=zb(l2(end))-drop; %new elevation of slope start
        l2=find(zb>=z2l & zb<=z1l); %update the hard layer
        
        
        peak=find(z==max(z(l2))); %find peak again
        
        %l2=find(zb>=z2l & zb<=z1l); %update the hard layer

        %place blocks every time hogback blocks fall
        blockstart=l2(end)+1; %reevaluate location of end of hogback/beginning of hillslope

        %choose block distribution pattern:
        %wdotdistnew=round(abs(stdev.*randn(1,nrocks))+blockstart); %use random distribution of block placement (right side of bell curve)
        %wdotdistnew=[blockstart:2:blockstart+((2*nrocks)-1)]; %place blocks immediately downslope, with one space in between each other
        wdotdistnew=[blockstart:blockstart+nrocks-1]; %place blocks immediately downslope, next to each other

        wdotdistnew(find(wdotdistnew>blockend))=[]; %make sure blocks don't go beyond calculation space
        zb(wdotdistnew)=zb(wdotdistnew)+blockheight+H(wdotdistnew); %add bedrock topography where there are blocks
        Hsavenew=H(wdotdistnew); %save soil height where there are now blocks
        Hsave=cat(2,Hsave,Hsavenew); %update saved soil height array

        H(wdotdistnew)=0; %set soil where block are to 0 (soil is saved, but must be 0 to allow new weathering rate)
        wdotdist=cat(2,wdotdist,wdotdistnew); %combine new block location array with old
        wdottracking=cat(2,wdottracking,zeros(size(wdotdistnew))); %add cells onto total weathering array to accomodate new blocks
        
    end


    %update weathering array to keep track of block heights
    if exponential==0
       
        %humped weathering rule
        wdot_Hsave=wdotnot.*exp((-(Hsave+(blockheight-wdottracking)))./wstar)+(((A*(Hsave+(blockheight-wdottracking)))./wstar).*exp(-(Hsave+(blockheight-wdottracking))./wstar));
        Hsave=Hsave+(wdot_Hsave*dt);
        %wdot(wdotdist)=wdotnot_hard.*exp(-H(wdotdist)./wstar)+(((A_l2*H(wdotdist))./wstar).*exp(-H(wdotdist))./wstar); %where there is a block, use hard wdot
        wdot(wdotdist) = wdotnot_hard;
        wdottracking=wdottracking+wdot(wdotdist).*dt; %keep track of total weathering
    
    end




    if exponential==1
       
        %exponential weathering rule
        wdot_Hsave=wdotnot.*exp(-(Hsave+(blockheight-wdottracking))./wstar);
        Hsave=Hsave+(wdot_Hsave*dt);
        wdot(wdotdist)=wdotnot_hard.*exp(-H(wdotdist)./wstar); %where there is a block, use hard wdot, allowing soil on top to affect weathering rate
        %wdot(wdotdist)=wdotnot_hard; %where there is a block, use hard wdot assuming no soil
        wdottracking=wdottracking+wdot(wdotdist).*dt; %keep track of total weathering
        
    
    end

    Hsave(find(wdotdist==length(x)-1))=[]; %remove saved soil height from array at old block locations 
    wdottracking(find(wdotdist==length(x)-1))=[]; %remove block weathering tracking from edge of domain
    wdotdist(find(wdotdist==length(x)-1))=[]; %remove block location tracking from edge of domain
    wdottracking(find(wdotdist==length(x)))=[]; %remove block weathering tracking from edge of domain
    wdotdist(find(wdotdist==length(x)))=[]; %remove block location tracking from edge of domain

    
    
    %Soil transport
    %calculate topographic slope
    slope=diff(z)/dx;

    % find peaks and troughs to find H at upslope cell
    signdiff = sign(slope);
    xedge = x(1:end-1)+(dx/2);
    lt = find(signdiff<0); % these will be slopes with peaks to their left
    rt = find(signdiff>0); % these will be slopes with peaks to their right
    
    %find soil thickness at upslope cell
    hedge = zeros(size(xedge));
    hedge(1) = H(2);
    hedge(end)= H(end-1);
    hedge(rt) = H(rt+1);
    hedge(lt) = H(lt);

    %linear flux from Johnstone and Hilley
    q = -k.*slope.*hstar.*(1-(exp(-hedge./hstar)));

    q = [q(1) q q(end)]; % boundary conditions
   
    
    % now conserve mass
    dh_dt=wdot-diff(q)/dx; %change in soil thickness per time step
    H=H+(dh_dt*dt); %new soil thickness 
    H = max(0,H); %prevent negative soil thickness
    
    if channel_on == 1
        zb(2:min(channel)-1) = zb(2:min(channel)-1)-(wdot(2:min(channel)-1));
    
    else
        zb(2:end)=zb(2:end)-(wdot(2:end)*dt); %update bedrock topography, except at boundaries
    end

    %Get rid of blocks once they have weathered to size cutoff
    H(wdotdist(find(wdottracking>=blockheight-cutoff)))=Hsave(find(wdottracking>=blockheight-cutoff));
    zb(wdotdist(find(wdottracking>=blockheight-cutoff)))=z(wdotdist(find(wdottracking>=blockheight-cutoff)))-H(wdotdist(find(wdottracking>=blockheight-cutoff)));
    Hsave(find(wdottracking>=blockheight-cutoff))=[];
    wdotdist(find(wdottracking>=blockheight-cutoff))=[];
    wdottracking(find(wdottracking>=blockheight-cutoff))=[];


    % now impose a channel boundary condition
    if(t(i)<=tchannel_stabilize)
       
        zb(1) = zb(1)-(edot_channel*dt);
        
        
        if channel_on ==1 %channel moves laterally with crest
            zb(channel) = zb(channel)-(edot_channel*dt);
            H(channel) = 0;
            
        else
            zb(end) = zb(end)-(edot_channel*dt);
        end
    
     else
        zb(1) = zb(1);
        zb(end) = zb(end);
        H(1)=1;
        H(end)=1;
    end
    
 
    z=H+zb; %update topography
    

    %Metrics
    if metric == 1
    
    height = blockheight-wdottracking; %current height of blocks
    wdotdistsort = sort(wdotdist); %locations of blocks
    wdotmean = wdot; %all weathering rates
    wdotmean(wdotdistsort) = NaN; %weathering rates of non-blocks only
    Hmean = H; %all soil depths
    Hmean(wdotdistsort) = NaN; %soil depth of non-blocks

        %calculate metrics after a certain amount of time
        if t(i)>1250000

            l = blockstart:1:blockstart+50;
            xmetrics = abs((x(l)) - (x(blockstart)));



            %metrics averaged over 4 cells
            for j = 1:length(l)


                  meanslope2(j) = ((z(l(j))-((sum(height(find(wdotdist == l(j)))))))-(z(l(j)+3)-(sum(height(find(wdotdist == l(j)+3))))))/(x(l(j)+3)-x(l(j)));
                  meanslope(j) = ((z(l(j))-(sum(height(find(wdotdist == l(j))))))-(z(l(j)+3)-(sum(height(find(wdotdist==l(j)+3))))))./((x(l(j)+3)-x(l(j))));
                  meanheight(j) = mean(height(find(wdotdist >= l(j) & wdotdist <= l(j)+3)));
                  %mean soil thickness between blocks
                  meanH(j) = nanmean(Hmean(l(j):l(j)+3));
                  meanwdot(j) = mean(wdot(l(j):l(j)+3));
                  %mean weathering rate in between blocks
                  meanwdot2(j) = nanmean(wdotmean(l(j):l(j)+3)); 
                  meanflux(j) = mean(q(l(j):l(j)+3));

            end

            meanheight(find(isnan(meanheight)==1)) = 0;


        end
    end



    %track sizes for fancy plotting
    markersize1 = wdotdist(find(wdottracking<=blockheight/4));
    markersize2 = wdotdist(find(wdottracking<=blockheight/2 & wdottracking>= blockheight/4));
    markersize3 = wdotdist(find(wdottracking<=3*blockheight/4 & wdottracking>= blockheight/2));
    markersize4 = wdotdist(find(wdottracking>=3*blockheight/4));

    if(rem(t(i),tplot) == 0)
        nframe = nframe+1;

    


    figure(1)
    subplot('position',[0.1 0.3 0.8 0.6]) 
    plot(x,z,'k')
    hold on
    plot(x,zb,'Color', [0.466667 0.533333 0.6])

    
    % now try to patch fill the hard stuff
    xhard = [x(1:peak) x(peak:l2(end)) fliplr(x(1:l2(end)))];
    zhard = [z1l(1:peak) zb(peak:l2(end)) fliplr(z2l(1:l2(end)))];
    fill(xhard,zhard, [0.698039 0.133333 0.133333]);

     plot(x(markersize1),z(markersize1),'sb','MarkerSize',8,'MarkerFaceColor',[0.698039 0.133333 0.133333],'MarkerEdgeColor','k')
     plot(x(markersize2),z(markersize2),'sb','MarkerSize',6,'MarkerFaceColor',[0.698039 0.133333 0.133333],'MarkerEdgeColor','k')
     plot(x(markersize3),z(markersize3),'sb','MarkerSize',3,'MarkerFaceColor',[0.698039 0.133333 0.133333],'MarkerEdgeColor','k')
     plot(x(markersize4),z(markersize4),'sb','MarkerSize',1,'MarkerFaceColor',[0.698039 0.133333 0.133333],'MarkerEdgeColor','k')

    set(gca,'fontsize',18,'fontname','arial');
    xlabel('Distance [m]','fontname','arial','fontsize',20) ;
    ylabel('Elevation [m]','fontname','arial','fontsize',20);
    axis([xlimfixed ylimfixed])
    ht=text(-230,max(ztopo0)+10,['  ',num2str(t(i)), ' years  '],'fontsize',21);
    hold off

    %  %subplot(2,1,2)
      subplot('position',[0.1 0.1 0.8 0.1])
      %plot(x,H,'-k')
      plot(x(wdotdist),blockheight-wdottracking,'.k','linewidth',2)
      xlabel('Distance [m]')
      ylabel('Block height [m]')
      axis([min(x) max(x)  0 blockheight])
    
    drawnow
            if (nframe==1)
                moviein(nplots,gcf);
            end
            M(:,nframe) = getframe(gcf);
 


     %store metrics to be averaged at end
     if metric == 1
         if t(i)>1250000


             meanslopesave=cat(1,meanslopesave,meanslope);
             meanslope2save = cat(1,meanslope2save,meanslope2);
             meanHsave=cat(1,meanHsave,meanH);
             meanheightsave=cat(1,meanheightsave,meanheight);
             meanwdotsave=cat(1,meanwdotsave,meanwdot);
             meanwdot2save=cat(1,meanwdot2save,meanwdot2);
             meanfluxsave=cat(1,meanfluxsave,meanflux);

         end
     end


    end
end




%movie2avi(M,'hogback_30deg_10m_conserve,expo','fps',8) % to convert M file to an avi movie
sorted=sort(wdotdist);
sorttrack=sort(wdottracking);

if metric == 1
    %analyticsave=[mean(analyticalslopesave(:,1)) mean(analyticalslopesave(:,2)) mean(analyticalslopesave(:,3)) mean(analyticalslopesave(:,4)) mean(analyticalslopesave(:,5)) mean(analyticalslopesave(:,6)) mean(analyticalslopesave(:,7)) mean(analyticalslopesave(:,8)) mean(analyticalslopesave(:,9)) mean(analyticalslopesave(:,10)) mean(analyticalslopesave(:,11)) mean(analyticalslopesave(:,12)) mean(analyticalslopesave(:,13)) mean(analyticalslopesave(:,14)) mean(analyticalslopesave(:,15)) mean(analyticalslopesave(:,16)) mean(analyticalslopesave(:,17)) mean(analyticalslopesave(:,18)) mean(analyticalslopesave(:,19)) mean(analyticalslopesave(:,20)) mean(analyticalslopesave(:,21))];
    modelsave=[mean(meanslopesave(:,1)) mean(meanslopesave(:,2)) mean(meanslopesave(:,3)) mean(meanslopesave(:,4)) mean(meanslopesave(:,5)) mean(meanslopesave(:,6)) mean(meanslopesave(:,7)) mean(meanslopesave(:,8)) mean(meanslopesave(:,9)) mean(meanslopesave(:,10)) mean(meanslopesave(:,11)) mean(meanslopesave(:,12)) mean(meanslopesave(:,13)) mean(meanslopesave(:,14)) mean(meanslopesave(:,15)) mean(meanslopesave(:,16)) mean(meanslopesave(:,17)) mean(meanslopesave(:,18)) mean(meanslopesave(:,19)) mean(meanslopesave(:,20)) mean(meanslopesave(:,21))];
    model2save = [mean(meanslope2save(:,1)) mean(meanslope2save(:,2)) mean(meanslope2save(:,3)) mean(meanslope2save(:,4)) mean(meanslope2save(:,5)) mean(meanslope2save(:,6)) mean(meanslope2save(:,7)) mean(meanslope2save(:,8)) mean(meanslope2save(:,9)) mean(meanslope2save(:,10)) mean(meanslope2save(:,11)) mean(meanslope2save(:,12)) mean(meanslope2save(:,13)) mean(meanslope2save(:,14)) mean(meanslope2save(:,15)) mean(meanslope2save(:,16)) mean(meanslope2save(:,17)) mean(meanslope2save(:,18)) mean(meanslope2save(:,19)) mean(meanslope2save(:,20)) mean(meanslope2save(:,21)) mean(meanslope2save(:,22)) mean(meanslope2save(:,23)) mean(meanslope2save(:,24)) mean(meanslope2save(:,25)) mean(meanslope2save(:,26)) mean(meanslope2save(:,27)) mean(meanslope2save(:,28)) mean(meanslope2save(:,29)) mean(meanslope2save(:,30)) mean(meanslope2save(:,31)) mean(meanslope2save(:,32)) mean(meanslope2save(:,33)) mean(meanslope2save(:,34)) mean(meanslope2save(:,35)) mean(meanslope2save(:,36)) mean(meanslope2save(:,37)) mean(meanslope2save(:,38)) mean(meanslope2save(:,39)) mean(meanslope2save(:,40)) mean(meanslope2save(:,41)) mean(meanslope2save(:,42)) mean(meanslope2save(:,43)) mean(meanslope2save(:,44)) mean(meanslope2save(:,45)) mean(meanslope2save(:,46)) mean(meanslope2save(:,47)) mean(meanslope2save(:,48)) mean(meanslope2save(:,49)) mean(meanslope2save(:,50))];
    heightsave=[mean(meanheightsave(:,1)) mean(meanheightsave(:,2)) mean(meanheightsave(:,3)) mean(meanheightsave(:,4)) mean(meanheightsave(:,5)) mean(meanheightsave(:,6)) mean(meanheightsave(:,7)) mean(meanheightsave(:,8)) mean(meanheightsave(:,9)) mean(meanheightsave(:,10)) mean(meanheightsave(:,11)) mean(meanheightsave(:,12)) mean(meanheightsave(:,13)) mean(meanheightsave(:,14)) mean(meanheightsave(:,15)) mean(meanheightsave(:,16)) mean(meanheightsave(:,17)) mean(meanheightsave(:,18)) mean(meanheightsave(:,19)) mean(meanheightsave(:,20)) mean(meanheightsave(:,21))];
    fluxsave=[mean(meanfluxsave(:,1)) mean(meanfluxsave(:,2)) mean(meanfluxsave(:,3)) mean(meanfluxsave(:,4)) mean(meanfluxsave(:,5)) mean(meanfluxsave(:,6)) mean(meanfluxsave(:,7)) mean(meanfluxsave(:,8)) mean(meanfluxsave(:,9)) mean(meanfluxsave(:,10)) mean(meanfluxsave(:,11)) mean(meanfluxsave(:,12)) mean(meanfluxsave(:,13)) mean(meanfluxsave(:,14)) mean(meanfluxsave(:,15)) mean(meanfluxsave(:,16)) mean(meanfluxsave(:,17)) mean(meanfluxsave(:,18)) mean(meanfluxsave(:,19)) mean(meanfluxsave(:,20)) mean(meanfluxsave(:,21))];
    Hsave=[mean(meanHsave(:,1)) mean(meanHsave(:,2)) mean(meanHsave(:,3)) mean(meanHsave(:,4)) mean(meanHsave(:,5)) mean(meanHsave(:,6)) mean(meanHsave(:,7)) mean(meanHsave(:,8)) mean(meanHsave(:,9)) mean(meanHsave(:,10)) mean(meanHsave(:,11)) mean(meanHsave(:,12)) mean(meanHsave(:,13)) mean(meanHsave(:,14)) mean(meanHsave(:,15)) mean(meanHsave(:,16)) mean(meanHsave(:,17)) mean(meanHsave(:,18)) mean(meanHsave(:,19)) mean(meanHsave(:,20)) mean(meanHsave(:,21))];
    wdotsave=[nanmean(meanwdotsave(:,1)) nanmean(meanwdotsave(:,2)) nanmean(meanwdotsave(:,3)) nanmean(meanwdotsave(:,4)) nanmean(meanwdotsave(:,5)) nanmean(meanwdotsave(:,6)) nanmean(meanwdotsave(:,7)) nanmean(meanwdotsave(:,8)) nanmean(meanwdotsave(:,9)) nanmean(meanwdotsave(:,10)) nanmean(meanwdotsave(:,11)) nanmean(meanwdotsave(:,12)) nanmean(meanwdotsave(:,13)) nanmean(meanwdotsave(:,14)) nanmean(meanwdotsave(:,15)) nanmean(meanwdotsave(:,16)) nanmean(meanwdotsave(:,17)) nanmean(meanwdotsave(:,18)) nanmean(meanwdotsave(:,19)) nanmean(meanwdotsave(:,20)) nanmean(meanwdotsave(:,21))];
    wdot2save=[nanmean(meanwdot2save(:,1)) nanmean(meanwdot2save(:,2)) nanmean(meanwdot2save(:,3)) nanmean(meanwdot2save(:,4)) nanmean(meanwdot2save(:,5)) nanmean(meanwdot2save(:,6)) nanmean(meanwdot2save(:,7)) nanmean(meanwdot2save(:,8)) nanmean(meanwdot2save(:,9)) nanmean(meanwdot2save(:,10)) nanmean(meanwdot2save(:,11)) nanmean(meanwdot2save(:,12)) nanmean(meanwdot2save(:,13)) nanmean(meanwdot2save(:,14)) nanmean(meanwdot2save(:,15)) nanmean(meanwdot2save(:,16)) nanmean(meanwdot2save(:,17)) nanmean(meanwdot2save(:,18)) nanmean(meanwdot2save(:,19)) nanmean(meanwdot2save(:,20)) nanmean(meanwdot2save(:,21))];



    figure
    plot(xmetrics(1:21),heightsave,'ok')
    xlabel('x')
    ylabel('block height [m]')

    figure
    plot(xmetrics(1:21),modelsave,'ok')
    xlabel('x')
    ylabel('model slope')

    figure
    plot(xmetrics(1:21),fluxsave,'ok')
    xlabel('x')
    ylabel('soil flux')

    figure
    plot(xmetrics(1:21),wdotsave,'ok')
    xlabel('x')
    ylabel('weathering rate')


    figure
    plot(xmetrics(1:21),wdot2save,'ok')
    xlabel('x')
    ylabel('weathering rate between blocks')


    figure
    plot(xmetrics(1:21),Hsave,'ok')
    xlabel('x')
    ylabel('soil thickness between blocks')
    
end




