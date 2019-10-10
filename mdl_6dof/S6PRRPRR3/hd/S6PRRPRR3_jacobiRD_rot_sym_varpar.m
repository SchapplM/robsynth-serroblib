% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:29
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRPRR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(6));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(6));
	t60 = cos(pkin(12));
	t58 = sin(pkin(12));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0, 0; 0, -t62 * t64, 0, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0, 0; 0, -t63 * t64, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:27
	% EndTime: 2019-10-09 22:29:27
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (78->41), mult. (285->88), div. (0->0), fcn. (301->10), ass. (0->35)
	t259 = sin(pkin(12));
	t262 = cos(pkin(12));
	t266 = sin(qJ(2));
	t264 = cos(pkin(6));
	t268 = cos(qJ(2));
	t283 = t264 * t268;
	t253 = -t259 * t266 + t262 * t283;
	t260 = sin(pkin(7));
	t261 = sin(pkin(6));
	t287 = t260 * t261;
	t263 = cos(pkin(7));
	t265 = sin(qJ(3));
	t286 = t263 * t265;
	t267 = cos(qJ(3));
	t285 = t263 * t267;
	t284 = t264 * t266;
	t282 = t265 * t266;
	t281 = t265 * t268;
	t280 = t266 * t267;
	t279 = t267 * t268;
	t277 = qJD(3) * t260 * t264;
	t276 = -t253 * t263 + t262 * t287;
	t274 = t259 * t283 + t262 * t266;
	t275 = -t259 * t287 + t263 * t274;
	t254 = t259 * t268 + t262 * t284;
	t273 = t259 * t284 - t262 * t268;
	t272 = -t263 * t279 + t282;
	t271 = -t263 * t280 - t281;
	t270 = -t263 * t281 - t280;
	t269 = t263 * t282 - t279;
	t252 = t273 * qJD(2);
	t251 = t274 * qJD(2);
	t250 = t254 * qJD(2);
	t249 = t253 * qJD(2);
	t1 = [0, t251 * t286 + t252 * t267 + (t265 * t274 + t273 * t285) * qJD(3), t252 * t285 + t251 * t265 + (t275 * t265 + t267 * t273) * qJD(3), 0, 0, 0; 0, -t249 * t286 - t250 * t267 + (-t253 * t265 - t254 * t285) * qJD(3), -t250 * t285 - t249 * t265 + (-t254 * t267 + t276 * t265) * qJD(3), 0, 0, 0; 0, (t270 * qJD(2) + t271 * qJD(3)) * t261, -t265 * t277 + (t271 * qJD(2) + t270 * qJD(3)) * t261, 0, 0, 0; 0, t251 * t285 - t252 * t265 + (t267 * t274 - t273 * t286) * qJD(3), -t252 * t286 + t251 * t267 + (-t265 * t273 + t275 * t267) * qJD(3), 0, 0, 0; 0, -t249 * t285 + t250 * t265 + (-t253 * t267 + t254 * t286) * qJD(3), t250 * t286 - t249 * t267 + (t254 * t265 + t276 * t267) * qJD(3), 0, 0, 0; 0, (t272 * qJD(2) + t269 * qJD(3)) * t261, -t267 * t277 + (t269 * qJD(2) + t272 * qJD(3)) * t261, 0, 0, 0; 0, -t251 * t260, 0, 0, 0, 0; 0, t249 * t260, 0, 0, 0, 0; 0, qJD(2) * t268 * t287, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:27
	% EndTime: 2019-10-09 22:29:28
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (132->51), mult. (453->101), div. (0->0), fcn. (493->12), ass. (0->38)
	t330 = sin(pkin(12));
	t334 = cos(pkin(12));
	t338 = sin(qJ(2));
	t336 = cos(pkin(6));
	t340 = cos(qJ(2));
	t348 = t336 * t340;
	t315 = -t330 * t338 + t334 * t348;
	t332 = sin(pkin(6));
	t352 = t330 * t332;
	t350 = t332 * t334;
	t349 = t336 * t338;
	t337 = sin(qJ(3));
	t347 = qJD(3) * t337;
	t339 = cos(qJ(3));
	t346 = qJD(3) * t339;
	t329 = sin(pkin(13));
	t333 = cos(pkin(13));
	t344 = t339 * t329 + t337 * t333;
	t321 = t337 * t329 - t339 * t333;
	t316 = t330 * t340 + t334 * t349;
	t343 = t330 * t348 + t334 * t338;
	t342 = t330 * t349 - t334 * t340;
	t341 = qJD(3) * t344;
	t319 = t321 * qJD(3);
	t335 = cos(pkin(7));
	t331 = sin(pkin(7));
	t320 = -t329 * t346 - t333 * t347;
	t314 = t342 * qJD(2);
	t313 = t343 * qJD(2);
	t312 = t316 * qJD(2);
	t311 = t315 * qJD(2);
	t310 = t344 * t335;
	t309 = t321 * t335;
	t308 = t335 * t341;
	t307 = (t329 * t347 - t333 * t346) * t335;
	t306 = t331 * t341;
	t305 = t331 * t319;
	t1 = [0, -t307 * t342 + t313 * t310 - t314 * t321 - t320 * t343, -t306 * t352 + t308 * t343 - t314 * t309 + t313 * t344 - t319 * t342, 0, 0, 0; 0, t316 * t307 - t311 * t310 + t312 * t321 + t315 * t320, t306 * t350 - t315 * t308 + t312 * t309 - t311 * t344 + t316 * t319, 0, 0, 0; 0, (t307 * t338 + t320 * t340 + (-t310 * t340 + t321 * t338) * qJD(2)) * t332, -t336 * t306 + (-t308 * t340 + t319 * t338 + (t309 * t338 - t340 * t344) * qJD(2)) * t332, 0, 0, 0; 0, -t308 * t342 - t313 * t309 - t314 * t344 - t319 * t343, t305 * t352 - t307 * t343 - t314 * t310 - t313 * t321 + t320 * t342, 0, 0, 0; 0, t316 * t308 + t311 * t309 + t312 * t344 + t315 * t319, -t305 * t350 + t315 * t307 + t312 * t310 + t311 * t321 - t316 * t320, 0, 0, 0; 0, (t308 * t338 + t319 * t340 + (t309 * t340 + t338 * t344) * qJD(2)) * t332, t336 * t305 + (t307 * t340 - t320 * t338 + (t310 * t338 + t321 * t340) * qJD(2)) * t332, 0, 0, 0; 0, -t313 * t331, 0, 0, 0, 0; 0, t311 * t331, 0, 0, 0, 0; 0, t332 * qJD(2) * t340 * t331, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:31
	% EndTime: 2019-10-09 22:29:32
	% DurationCPUTime: 0.97s
	% Computational Cost: add. (428->105), mult. (1395->201), div. (0->0), fcn. (1591->14), ass. (0->73)
	t550 = sin(pkin(13));
	t558 = sin(qJ(3));
	t591 = cos(pkin(13));
	t592 = cos(qJ(3));
	t563 = -t558 * t550 + t592 * t591;
	t595 = t563 * qJD(3);
	t542 = -t592 * t550 - t558 * t591;
	t594 = qJD(3) * t542;
	t551 = sin(pkin(12));
	t554 = cos(pkin(12));
	t559 = sin(qJ(2));
	t556 = cos(pkin(6));
	t561 = cos(qJ(2));
	t582 = t556 * t561;
	t535 = -t551 * t559 + t554 * t582;
	t553 = sin(pkin(6));
	t590 = t551 * t553;
	t552 = sin(pkin(7));
	t588 = t552 * t553;
	t557 = sin(qJ(5));
	t587 = t552 * t557;
	t560 = cos(qJ(5));
	t586 = t552 * t560;
	t585 = t553 * t554;
	t555 = cos(pkin(7));
	t584 = t553 * t555;
	t583 = t556 * t559;
	t580 = qJD(5) * t557;
	t579 = qJD(5) * t560;
	t578 = t559 * t588;
	t575 = qJD(2) * t588;
	t572 = t559 * t575;
	t571 = t561 * t575;
	t528 = t563 * t555;
	t569 = t528 * t561 + t542 * t559;
	t529 = t542 * t555;
	t568 = -t529 * t561 + t559 * t563;
	t567 = t529 * t559 + t561 * t563;
	t536 = t551 * t561 + t554 * t583;
	t565 = t551 * t582 + t554 * t559;
	t564 = t551 * t583 - t554 * t561;
	t522 = t595 * t552;
	t524 = t595 * t555;
	t530 = t535 * qJD(2);
	t531 = t536 * qJD(2);
	t503 = -t522 * t585 + t535 * t524 + t531 * t529 + t530 * t563 + t536 * t594;
	t532 = t565 * qJD(2);
	t533 = t564 * qJD(2);
	t505 = t522 * t590 - t524 * t565 - t533 * t529 - t532 * t563 - t564 * t594;
	t509 = t556 * t522 + (t567 * qJD(2) + t524 * t561 + t559 * t594) * t553;
	t534 = t556 * t555 - t561 * t588;
	t527 = t542 * t552;
	t526 = t563 * t552;
	t525 = t555 * t594;
	t523 = t552 * t594;
	t521 = t551 * t584 + t552 * t565;
	t520 = -t535 * t552 - t554 * t584;
	t519 = t567 * t553;
	t518 = -t529 * t564 - t563 * t565;
	t517 = t536 * t529 + t535 * t563;
	t516 = -t556 * t527 + t568 * t553;
	t515 = t556 * t526 + t569 * t553;
	t514 = -t527 * t590 + t529 * t565 - t563 * t564;
	t513 = t526 * t590 - t528 * t565 - t542 * t564;
	t512 = t527 * t585 - t535 * t529 + t536 * t563;
	t511 = -t526 * t585 + t535 * t528 + t536 * t542;
	t510 = (-t568 * qJD(2) - t524 * t559 + t561 * t594) * t553;
	t508 = t556 * t523 + (t525 * t561 - t595 * t559 + (-t528 * t559 + t542 * t561) * qJD(2)) * t553;
	t507 = t524 * t564 - t532 * t529 + t533 * t563 - t565 * t594;
	t506 = -t536 * t524 + t530 * t529 - t531 * t563 + t535 * t594;
	t504 = t523 * t590 - t525 * t565 + t533 * t528 - t532 * t542 + t564 * t595;
	t502 = -t523 * t585 + t535 * t525 - t531 * t528 + t530 * t542 - t536 * t595;
	t1 = [0, -t532 * t587 + t507 * t560 + (-t518 * t557 - t564 * t586) * qJD(5), t504 * t560 - t513 * t580, 0, -t533 * t586 - t505 * t557 + (-t514 * t560 - t521 * t557) * qJD(5), 0; 0, t530 * t587 + t506 * t560 + (-t517 * t557 + t536 * t586) * qJD(5), t502 * t560 - t511 * t580, 0, t531 * t586 - t503 * t557 + (-t512 * t560 - t520 * t557) * qJD(5), 0; 0, t557 * t571 + t510 * t560 + (-t519 * t557 + t560 * t578) * qJD(5), t508 * t560 - t515 * t580, 0, t560 * t572 - t509 * t557 + (-t516 * t560 - t534 * t557) * qJD(5), 0; 0, -t532 * t586 - t507 * t557 + (-t518 * t560 + t564 * t587) * qJD(5), -t504 * t557 - t513 * t579, 0, t533 * t587 - t505 * t560 + (t514 * t557 - t521 * t560) * qJD(5), 0; 0, t530 * t586 - t506 * t557 + (-t517 * t560 - t536 * t587) * qJD(5), -t502 * t557 - t511 * t579, 0, -t531 * t587 - t503 * t560 + (t512 * t557 - t520 * t560) * qJD(5), 0; 0, t560 * t571 - t510 * t557 + (-t519 * t560 - t557 * t578) * qJD(5), -t508 * t557 - t515 * t579, 0, -t557 * t572 - t509 * t560 + (t516 * t557 - t534 * t560) * qJD(5), 0; 0, -t525 * t564 - t532 * t528 - t533 * t542 - t565 * t595, t505, 0, 0, 0; 0, t536 * t525 + t530 * t528 + t531 * t542 + t535 * t595, t503, 0, 0, 0; 0, (t569 * qJD(2) + t525 * t559 + t561 * t595) * t553, t509, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:35
	% EndTime: 2019-10-09 22:29:36
	% DurationCPUTime: 1.80s
	% Computational Cost: add. (1228->163), mult. (3848->302), div. (0->0), fcn. (4524->16), ass. (0->111)
	t747 = sin(pkin(13));
	t756 = sin(qJ(3));
	t806 = cos(pkin(13));
	t807 = cos(qJ(3));
	t771 = -t756 * t747 + t807 * t806;
	t810 = t771 * qJD(3);
	t739 = -t807 * t747 - t756 * t806;
	t809 = qJD(3) * t739;
	t748 = sin(pkin(12));
	t750 = sin(pkin(6));
	t805 = t748 * t750;
	t749 = sin(pkin(7));
	t804 = t749 * t750;
	t755 = sin(qJ(5));
	t803 = t749 * t755;
	t759 = cos(qJ(5));
	t802 = t749 * t759;
	t751 = cos(pkin(12));
	t801 = t750 * t751;
	t752 = cos(pkin(7));
	t800 = t750 * t752;
	t753 = cos(pkin(6));
	t757 = sin(qJ(2));
	t799 = t753 * t757;
	t760 = cos(qJ(2));
	t798 = t753 * t760;
	t796 = qJD(2) * t757;
	t795 = qJD(5) * t755;
	t794 = qJD(5) * t759;
	t754 = sin(qJ(6));
	t793 = qJD(6) * t754;
	t758 = cos(qJ(6));
	t792 = qJD(6) * t758;
	t791 = qJD(6) * t759;
	t790 = t757 * t804;
	t789 = t760 * t804;
	t788 = t751 * t798;
	t786 = t750 * t796;
	t783 = t749 * t786;
	t782 = qJD(2) * t789;
	t723 = t771 * t749;
	t725 = t771 * t752;
	t732 = -t748 * t757 + t788;
	t733 = t748 * t760 + t751 * t799;
	t697 = -t723 * t801 + t732 * t725 + t733 * t739;
	t719 = t810 * t749;
	t721 = t810 * t752;
	t726 = t739 * t752;
	t727 = -qJD(2) * t788 + t748 * t796;
	t728 = t733 * qJD(2);
	t763 = -t719 * t801 + t732 * t721 + t728 * t726 - t727 * t771 + t733 * t809;
	t780 = -t697 * t791 + t763;
	t772 = t748 * t799 - t751 * t760;
	t773 = t748 * t798 + t751 * t757;
	t700 = t723 * t805 - t725 * t773 - t739 * t772;
	t729 = t773 * qJD(2);
	t730 = t772 * qJD(2);
	t762 = t719 * t805 - t721 * t773 - t730 * t726 - t729 * t771 - t772 * t809;
	t779 = -t700 * t791 + t762;
	t777 = t725 * t760 + t739 * t757;
	t705 = t753 * t723 + t777 * t750;
	t775 = t726 * t757 + t760 * t771;
	t761 = t753 * t719 + (t775 * qJD(2) + t721 * t760 + t757 * t809) * t750;
	t778 = -t705 * t791 + t761;
	t714 = -t732 * t749 - t751 * t800;
	t724 = t739 * t749;
	t770 = t724 * t801 - t732 * t726 + t733 * t771;
	t686 = t714 * t755 + t759 * t770;
	t685 = t714 * t759 - t755 * t770;
	t715 = t748 * t800 + t749 * t773;
	t769 = -t724 * t805 + t726 * t773 - t771 * t772;
	t688 = t715 * t755 + t759 * t769;
	t687 = t715 * t759 - t755 * t769;
	t731 = t753 * t752 - t789;
	t776 = -t726 * t760 + t757 * t771;
	t764 = -t753 * t724 + t776 * t750;
	t703 = t731 * t755 + t759 * t764;
	t702 = t731 * t759 - t755 * t764;
	t708 = t733 * t726 + t732 * t771;
	t694 = t708 * t759 + t733 * t803;
	t710 = -t726 * t772 - t771 * t773;
	t695 = t710 * t759 - t772 * t803;
	t713 = t775 * t750;
	t711 = t713 * t759 + t755 * t790;
	t720 = t749 * t809;
	t722 = t752 * t809;
	t676 = -t720 * t801 + t732 * t722 - t728 * t725 - t727 * t739 - t733 * t810;
	t768 = -qJD(6) * t770 - t676 * t759 + t697 * t795;
	t679 = t720 * t805 - t722 * t773 + t730 * t725 - t729 * t739 + t772 * t810;
	t767 = -qJD(6) * t769 - t679 * t759 + t700 * t795;
	t690 = t753 * t720 - t725 * t786 + (-t810 * t757 + (qJD(2) * t739 + t722) * t760) * t750;
	t766 = -qJD(6) * t764 - t690 * t759 + t705 * t795;
	t712 = (t725 * t757 - t739 * t760) * t750;
	t709 = -t725 * t772 + t739 * t773;
	t707 = t733 * t725 - t732 * t739;
	t693 = (-t776 * qJD(2) - t721 * t757 + t760 * t809) * t750;
	t692 = (t777 * qJD(2) + t722 * t757 + t760 * t810) * t750;
	t684 = t721 * t772 - t729 * t726 + t730 * t771 - t773 * t809;
	t683 = -t722 * t772 - t729 * t725 - t730 * t739 - t773 * t810;
	t682 = -t733 * t721 - t727 * t726 - t728 * t771 + t732 * t809;
	t681 = t733 * t722 - t727 * t725 + t728 * t739 + t732 * t810;
	t674 = t755 * t782 + t693 * t759 + (-t713 * t755 + t759 * t790) * qJD(5);
	t673 = t702 * qJD(5) + t755 * t783 + t759 * t761;
	t672 = -t703 * qJD(5) - t755 * t761 + t759 * t783;
	t671 = -t729 * t803 + t684 * t759 + (-t710 * t755 - t772 * t802) * qJD(5);
	t670 = -t727 * t803 + t682 * t759 + (-t708 * t755 + t733 * t802) * qJD(5);
	t669 = t687 * qJD(5) - t730 * t803 + t759 * t762;
	t668 = -t688 * qJD(5) - t730 * t802 - t755 * t762;
	t667 = t685 * qJD(5) + t728 * t803 + t759 * t763;
	t666 = -t686 * qJD(5) + t728 * t802 - t755 * t763;
	t1 = [0, t671 * t758 + t683 * t754 + (-t695 * t754 + t709 * t758) * qJD(6), t779 * t754 - t767 * t758, 0, t668 * t758 - t687 * t793, -t669 * t754 - t679 * t758 + (-t688 * t758 + t700 * t754) * qJD(6); 0, t670 * t758 + t681 * t754 + (-t694 * t754 + t707 * t758) * qJD(6), t780 * t754 - t768 * t758, 0, t666 * t758 - t685 * t793, -t667 * t754 - t676 * t758 + (-t686 * t758 + t697 * t754) * qJD(6); 0, t674 * t758 + t692 * t754 + (-t711 * t754 + t712 * t758) * qJD(6), t778 * t754 - t766 * t758, 0, t672 * t758 - t702 * t793, -t673 * t754 - t690 * t758 + (-t703 * t758 + t705 * t754) * qJD(6); 0, -t671 * t754 + t683 * t758 + (-t695 * t758 - t709 * t754) * qJD(6), t767 * t754 + t779 * t758, 0, -t668 * t754 - t687 * t792, -t669 * t758 + t679 * t754 + (t688 * t754 + t700 * t758) * qJD(6); 0, -t670 * t754 + t681 * t758 + (-t694 * t758 - t707 * t754) * qJD(6), t768 * t754 + t780 * t758, 0, -t666 * t754 - t685 * t792, -t667 * t758 + t676 * t754 + (t686 * t754 + t697 * t758) * qJD(6); 0, -t674 * t754 + t692 * t758 + (-t711 * t758 - t712 * t754) * qJD(6), t766 * t754 + t778 * t758, 0, -t672 * t754 - t702 * t792, -t673 * t758 + t690 * t754 + (t703 * t754 + t705 * t758) * qJD(6); 0, t695 * qJD(5) + t684 * t755 + t729 * t802, t679 * t755 + t700 * t794, 0, t669, 0; 0, t694 * qJD(5) + t682 * t755 + t727 * t802, t676 * t755 + t697 * t794, 0, t667, 0; 0, t711 * qJD(5) + t693 * t755 - t759 * t782, t690 * t755 + t705 * t794, 0, t673, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end