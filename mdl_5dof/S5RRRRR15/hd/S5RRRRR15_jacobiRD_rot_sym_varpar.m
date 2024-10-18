% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: S5RRRRR15_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRRR15_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR15_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (27->13), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t158 = sin(qJ(2));
	t159 = sin(qJ(1));
	t171 = t158 * t159;
	t161 = cos(qJ(1));
	t170 = t158 * t161;
	t160 = cos(qJ(2));
	t169 = t159 * t160;
	t168 = t160 * t161;
	t156 = sin(pkin(5));
	t167 = qJD(1) * t156;
	t166 = qJD(2) * t156;
	t157 = cos(pkin(5));
	t165 = -t157 * t168 + t171;
	t164 = t157 * t169 + t170;
	t163 = t157 * t170 + t169;
	t162 = t157 * t171 - t168;
	t155 = t162 * qJD(1) + t165 * qJD(2);
	t154 = t164 * qJD(1) + t163 * qJD(2);
	t153 = t163 * qJD(1) + t164 * qJD(2);
	t152 = t165 * qJD(1) + t162 * qJD(2);
	t1 = [t155, t152, 0, 0, 0; -t153, -t154, 0, 0, 0; 0, -t158 * t166, 0, 0, 0; t154, t153, 0, 0, 0; t152, t155, 0, 0, 0; 0, -t160 * t166, 0, 0, 0; -t159 * t167, 0, 0, 0, 0; t161 * t167, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:08
	% EndTime: 2024-09-27 22:28:08
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (245->23), mult. (152->31), div. (0->0), fcn. (132->9), ass. (0->31)
	t282 = sin(qJ(1));
	t295 = t282 / 0.2e1;
	t280 = qJ(2) + qJ(3);
	t275 = pkin(5) + t280;
	t271 = sin(t275);
	t276 = pkin(5) - t280;
	t272 = sin(t276);
	t268 = t271 - t272;
	t278 = cos(t280);
	t283 = cos(qJ(1));
	t285 = t282 * t278 + t283 * t268 / 0.2e1;
	t296 = t268 * t295 - t283 * t278;
	t274 = cos(t276);
	t279 = qJD(2) + qJD(3);
	t270 = t279 * t274;
	t273 = cos(t275);
	t293 = t279 * t273;
	t266 = t270 + t293;
	t277 = sin(t280);
	t289 = t283 * t277;
	t260 = t285 * qJD(1) + t266 * t295 + t279 * t289;
	t294 = -t283 / 0.2e1;
	t292 = t282 * t277;
	t287 = t274 + t273;
	t286 = qJD(1) * sin(pkin(5));
	t262 = t296 * qJD(1) + t266 * t294 + t279 * t292;
	t265 = t293 / 0.2e1 - t270 / 0.2e1;
	t264 = (-t271 / 0.2e1 - t272 / 0.2e1) * t279;
	t261 = (t287 * t295 + t289) * qJD(1) + t285 * t279;
	t259 = (t287 * t294 + t292) * qJD(1) + t296 * t279;
	t1 = [t262, t259, t259, 0, 0; -t260, -t261, -t261, 0, 0; 0, t265, t265, 0, 0; t261, t260, t260, 0, 0; t259, t262, t262, 0, 0; 0, t264, t264, 0, 0; -t282 * t286, 0, 0, 0, 0; t283 * t286, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:08
	% EndTime: 2024-09-27 22:28:08
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (499->25), mult. (236->34), div. (0->0), fcn. (176->9), ass. (0->29)
	t297 = qJ(2) + qJ(3) + qJ(4);
	t293 = pkin(5) - t297;
	t312 = cos(t293);
	t292 = pkin(5) + t297;
	t313 = cos(t292) / 0.2e1;
	t285 = t312 / 0.2e1 + t313;
	t296 = qJD(2) + qJD(3) + qJD(4);
	t283 = t285 * t296;
	t295 = cos(t297);
	t299 = sin(qJ(1));
	t300 = cos(qJ(1));
	t305 = -sin(t292) / 0.2e1;
	t311 = sin(t293);
	t303 = t311 / 0.2e1 + t305;
	t294 = sin(t297);
	t310 = t294 * t300;
	t277 = (t295 * t299 - t300 * t303) * qJD(1) + t299 * t283 + t296 * t310;
	t309 = t296 * t299;
	t308 = t300 * t295;
	t307 = qJD(1) * t299;
	t306 = qJD(1) * t300;
	t301 = t303 * t296;
	t279 = -t300 * t283 + t294 * t309 + (-t299 * t303 - t308) * qJD(1);
	t298 = sin(pkin(5));
	t282 = (t313 - t312 / 0.2e1) * t296;
	t281 = (t305 - t311 / 0.2e1) * t296;
	t278 = t295 * t309 - t300 * t301 + (t285 * t299 + t310) * qJD(1);
	t276 = -t285 * t306 + t294 * t307 - t296 * t308 - t299 * t301;
	t1 = [t279, t276, t276, t276, 0; -t277, -t278, -t278, -t278, 0; 0, t282, t282, t282, 0; t278, t277, t277, t277, 0; t276, t279, t279, t279, 0; 0, t281, t281, t281, 0; -t298 * t307, 0, 0, 0, 0; t298 * t306, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:09
	% EndTime: 2024-09-27 22:28:09
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (732->64), mult. (1069->124), div. (0->0), fcn. (1109->10), ass. (0->69)
	t591 = qJ(2) + qJ(3) + qJ(4);
	t588 = sin(t591);
	t589 = cos(t591);
	t596 = sin(qJ(5));
	t598 = cos(qJ(5));
	t594 = cos(pkin(6));
	t595 = cos(pkin(5));
	t597 = sin(qJ(1));
	t632 = t597 * t596;
	t623 = t595 * t632;
	t599 = cos(qJ(1));
	t629 = t599 * t598;
	t578 = t594 * t623 - t629;
	t621 = t595 * t629;
	t582 = t594 * t632 - t621;
	t585 = -t594 * t629 + t623;
	t590 = qJD(2) + qJD(3) + qJD(4);
	t613 = t582 * qJD(1) + t585 * qJD(5) + t578 * t590;
	t630 = t599 * t596;
	t624 = t595 * t630;
	t631 = t598 * t597;
	t580 = t594 * t624 + t631;
	t622 = t595 * t631;
	t584 = t594 * t630 + t622;
	t602 = t594 * t622 + t630;
	t615 = t580 * qJD(1) + t602 * qJD(5) + t584 * t590;
	t626 = qJD(5) * t597;
	t627 = qJD(1) * t599;
	t592 = sin(pkin(6));
	t593 = sin(pkin(5));
	t637 = t592 * t593;
	t639 = t613 * t588 - t615 * t589 + (t596 * t627 + t598 * t626) * t637;
	t581 = t594 * t621 - t632;
	t614 = -t585 * qJD(1) - t582 * qJD(5) + t581 * t590;
	t583 = t594 * t631 + t624;
	t616 = t602 * qJD(1) + t580 * qJD(5) + t583 * t590;
	t625 = qJD(5) * t599;
	t628 = qJD(1) * t597;
	t638 = t614 * t588 + t616 * t589 - (t596 * t625 + t598 * t628) * t637;
	t636 = t594 * t596;
	t635 = t594 * t598;
	t634 = t595 * t597;
	t633 = t595 * t599;
	t620 = qJD(1) * t593 * t594;
	t619 = qJD(5) * t592 * t595;
	t618 = -t581 * qJD(1) + t578 * qJD(5) + t585 * t590;
	t617 = t578 * qJD(1) - t581 * qJD(5) + t582 * t590;
	t612 = t583 * qJD(1) + t584 * qJD(5) + t590 * t602;
	t611 = t584 * qJD(1) + t583 * qJD(5) + t580 * t590;
	t610 = t588 * t596 - t589 * t635;
	t609 = -t588 * t598 - t589 * t636;
	t608 = -t588 * t635 - t589 * t596;
	t607 = t588 * t636 - t589 * t598;
	t606 = -t588 * t597 + t589 * t633;
	t605 = -t588 * t599 - t589 * t634;
	t604 = -t588 * t633 - t589 * t597;
	t603 = -t588 * t634 + t589 * t599;
	t586 = t590 * t589 * t637;
	t569 = (t608 * qJD(5) + t609 * t590) * t593;
	t568 = (t607 * qJD(5) + t610 * t590) * t593;
	t567 = (t603 * qJD(1) + t606 * t590) * t592;
	t566 = (t604 * qJD(1) + t605 * t590) * t592;
	t565 = t615 * t588 + t613 * t589;
	t564 = t616 * t588 - t614 * t589;
	t563 = t617 * t588 - t611 * t589;
	t562 = -t618 * t588 + t612 * t589;
	t561 = (-t596 * t628 + t598 * t625) * t637 + t617 * t589 + t611 * t588;
	t560 = (-t596 * t626 + t598 * t627) * t637 + t618 * t589 + t612 * t588;
	t1 = [t561, t565, t565, t565, t560; t639, t563, t563, t563, -t638; 0, t569, t569, t569, -t596 * t619 + (t609 * qJD(5) + t608 * t590) * t593; t638, t562, t562, t562, -t639; t560, t564, t564, t564, t561; 0, t568, t568, t568, -t598 * t619 + (t610 * qJD(5) + t607 * t590) * t593; -t597 * t620 + (t605 * qJD(1) + t604 * t590) * t592, t566, t566, t566, 0; t599 * t620 + (t606 * qJD(1) + t603 * t590) * t592, t567, t567, t567, 0; 0, t586, t586, t586, 0;];
	JRD_rot = t1;
end