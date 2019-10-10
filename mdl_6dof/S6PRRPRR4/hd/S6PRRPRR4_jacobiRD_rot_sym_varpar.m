% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:31
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRPRR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(6));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(6));
	t60 = cos(pkin(11));
	t58 = sin(pkin(11));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0, 0; 0, -t62 * t64, 0, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0, 0; 0, -t63 * t64, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:30
	% EndTime: 2019-10-09 22:31:30
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (36->22), mult. (140->59), div. (0->0), fcn. (148->8), ass. (0->26)
	t211 = sin(pkin(6));
	t214 = sin(qJ(3));
	t227 = t211 * t214;
	t216 = cos(qJ(3));
	t226 = t211 * t216;
	t213 = cos(pkin(6));
	t215 = sin(qJ(2));
	t225 = t213 * t215;
	t217 = cos(qJ(2));
	t224 = t213 * t217;
	t223 = qJD(2) * t215;
	t222 = qJD(3) * t214;
	t221 = qJD(3) * t216;
	t220 = qJD(3) * t217;
	t219 = t211 * qJD(2) * t217;
	t210 = sin(pkin(11));
	t212 = cos(pkin(11));
	t206 = -t210 * t215 + t212 * t224;
	t207 = t210 * t217 + t212 * t225;
	t208 = -t210 * t224 - t212 * t215;
	t218 = t210 * t225 - t212 * t217;
	t205 = t218 * qJD(2);
	t204 = t208 * qJD(2);
	t203 = t207 * qJD(2);
	t202 = t206 * qJD(2);
	t1 = [0, t205 * t216 - t208 * t222, -t204 * t214 + (-t210 * t227 + t216 * t218) * qJD(3), 0, 0, 0; 0, -t203 * t216 - t206 * t222, -t202 * t214 + (-t207 * t216 + t212 * t227) * qJD(3), 0, 0, 0; 0, (-t214 * t220 - t216 * t223) * t211, -t214 * t219 + (-t213 * t214 - t215 * t226) * qJD(3), 0, 0, 0; 0, -t205 * t214 - t208 * t221, -t204 * t216 + (-t210 * t226 - t214 * t218) * qJD(3), 0, 0, 0; 0, t203 * t214 - t206 * t221, -t202 * t216 + (t207 * t214 + t212 * t226) * qJD(3), 0, 0, 0; 0, (t214 * t223 - t216 * t220) * t211, -t216 * t219 + (-t213 * t216 + t215 * t227) * qJD(3), 0, 0, 0; 0, t204, 0, 0, 0, 0; 0, t202, 0, 0, 0, 0; 0, t219, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:30
	% EndTime: 2019-10-09 22:31:31
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (36->22), mult. (140->59), div. (0->0), fcn. (148->8), ass. (0->26)
	t269 = sin(pkin(6));
	t272 = sin(qJ(3));
	t285 = t269 * t272;
	t274 = cos(qJ(3));
	t284 = t269 * t274;
	t271 = cos(pkin(6));
	t273 = sin(qJ(2));
	t283 = t271 * t273;
	t275 = cos(qJ(2));
	t282 = t271 * t275;
	t281 = qJD(2) * t273;
	t280 = qJD(3) * t272;
	t279 = qJD(3) * t274;
	t278 = qJD(3) * t275;
	t277 = t269 * qJD(2) * t275;
	t268 = sin(pkin(11));
	t270 = cos(pkin(11));
	t264 = -t268 * t273 + t270 * t282;
	t265 = t268 * t275 + t270 * t283;
	t266 = -t268 * t282 - t270 * t273;
	t276 = t268 * t283 - t270 * t275;
	t263 = t276 * qJD(2);
	t262 = t266 * qJD(2);
	t261 = t265 * qJD(2);
	t260 = t264 * qJD(2);
	t1 = [0, t263 * t274 - t266 * t280, -t262 * t272 + (-t268 * t285 + t274 * t276) * qJD(3), 0, 0, 0; 0, -t261 * t274 - t264 * t280, -t260 * t272 + (-t265 * t274 + t270 * t285) * qJD(3), 0, 0, 0; 0, (-t272 * t278 - t274 * t281) * t269, -t272 * t277 + (-t271 * t272 - t273 * t284) * qJD(3), 0, 0, 0; 0, t262, 0, 0, 0, 0; 0, t260, 0, 0, 0, 0; 0, t277, 0, 0, 0, 0; 0, t263 * t272 + t266 * t279, t262 * t274 + (t268 * t284 + t272 * t276) * qJD(3), 0, 0, 0; 0, -t261 * t272 + t264 * t279, t260 * t274 + (-t265 * t272 - t270 * t284) * qJD(3), 0, 0, 0; 0, (-t272 * t281 + t274 * t278) * t269, t274 * t277 + (t271 * t274 - t273 * t285) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:31
	% EndTime: 2019-10-09 22:31:31
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (213->51), mult. (692->93), div. (0->0), fcn. (764->10), ass. (0->49)
	t324 = sin(pkin(11));
	t326 = cos(pkin(11));
	t333 = cos(qJ(2));
	t327 = cos(pkin(6));
	t330 = sin(qJ(2));
	t358 = t327 * t330;
	t315 = t324 * t333 + t326 * t358;
	t332 = cos(qJ(3));
	t325 = sin(pkin(6));
	t329 = sin(qJ(3));
	t360 = t325 * t329;
	t305 = t315 * t332 - t326 * t360;
	t368 = qJD(3) - qJD(5);
	t359 = t325 * t332;
	t319 = t327 * t329 + t330 * t359;
	t353 = t325 * qJD(2) * t333;
	t308 = t319 * qJD(3) + t329 * t353;
	t318 = -t327 * t332 + t330 * t360;
	t309 = -t318 * qJD(3) + t332 * t353;
	t328 = sin(qJ(5));
	t331 = cos(qJ(5));
	t367 = -t308 * t331 + t309 * t328 + (t318 * t328 + t319 * t331) * qJD(5);
	t366 = t308 * t328 + t309 * t331 + (t318 * t331 - t319 * t328) * qJD(5);
	t336 = t324 * t358 - t326 * t333;
	t307 = t324 * t360 - t332 * t336;
	t357 = t327 * t333;
	t337 = t324 * t357 + t326 * t330;
	t312 = t337 * qJD(2);
	t302 = t307 * qJD(3) - t312 * t329;
	t338 = t324 * t359 + t329 * t336;
	t303 = t338 * qJD(3) - t312 * t332;
	t365 = (t307 * t331 - t328 * t338) * qJD(5) - t302 * t331 + t303 * t328;
	t364 = (-t307 * t328 - t331 * t338) * qJD(5) + t302 * t328 + t303 * t331;
	t354 = t326 * t357;
	t356 = qJD(2) * t330;
	t310 = -qJD(2) * t354 + t324 * t356;
	t300 = t305 * qJD(3) - t310 * t329;
	t304 = t315 * t329 + t326 * t359;
	t301 = -t304 * qJD(3) - t310 * t332;
	t363 = -t300 * t331 + t301 * t328 + (t304 * t328 + t305 * t331) * qJD(5);
	t362 = t300 * t328 + t301 * t331 + (t304 * t331 - t305 * t328) * qJD(5);
	t340 = t328 * t332 - t329 * t331;
	t339 = t328 * t329 + t331 * t332;
	t335 = t368 * t340;
	t334 = t368 * t339;
	t314 = -t324 * t330 + t354;
	t313 = t336 * qJD(2);
	t311 = t315 * qJD(2);
	t1 = [0, t339 * t313 - t335 * t337, t365, 0, -t365, 0; 0, -t339 * t311 + t335 * t314, t363, 0, -t363, 0; 0, (t335 * t333 - t339 * t356) * t325, t367, 0, -t367, 0; 0, -t340 * t313 - t334 * t337, t364, 0, -t364, 0; 0, t340 * t311 + t334 * t314, t362, 0, -t362, 0; 0, (t334 * t333 + t340 * t356) * t325, t366, 0, -t366, 0; 0, t312, 0, 0, 0, 0; 0, t310, 0, 0, 0, 0; 0, -t353, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:35
	% EndTime: 2019-10-09 22:31:36
	% DurationCPUTime: 0.84s
	% Computational Cost: add. (549->95), mult. (1712->173), div. (0->0), fcn. (1948->12), ass. (0->74)
	t624 = sin(pkin(11));
	t626 = cos(pkin(11));
	t635 = cos(qJ(2));
	t627 = cos(pkin(6));
	t631 = sin(qJ(2));
	t657 = t627 * t631;
	t614 = t624 * t635 + t626 * t657;
	t630 = sin(qJ(3));
	t625 = sin(pkin(6));
	t634 = cos(qJ(3));
	t659 = t625 * t634;
	t602 = t614 * t630 + t626 * t659;
	t661 = t625 * t630;
	t603 = t614 * t634 - t626 * t661;
	t629 = sin(qJ(5));
	t633 = cos(qJ(5));
	t589 = t602 * t629 + t603 * t633;
	t656 = t627 * t635;
	t651 = t626 * t656;
	t655 = qJD(2) * t631;
	t609 = -qJD(2) * t651 + t624 * t655;
	t593 = t603 * qJD(3) - t609 * t630;
	t594 = -t602 * qJD(3) - t609 * t634;
	t575 = t589 * qJD(5) - t593 * t633 + t594 * t629;
	t588 = t602 * t633 - t603 * t629;
	t628 = sin(qJ(6));
	t632 = cos(qJ(6));
	t653 = qJD(6) * t632;
	t675 = t575 * t628 - t588 * t653;
	t654 = qJD(6) * t628;
	t674 = t575 * t632 + t588 * t654;
	t638 = t624 * t657 - t626 * t635;
	t605 = t624 * t661 - t634 * t638;
	t640 = t624 * t659 + t630 * t638;
	t592 = t605 * t633 - t629 * t640;
	t639 = t624 * t656 + t626 * t631;
	t611 = t639 * qJD(2);
	t595 = t605 * qJD(3) - t611 * t630;
	t596 = t640 * qJD(3) - t611 * t634;
	t578 = t592 * qJD(5) - t595 * t633 + t596 * t629;
	t591 = -t605 * t629 - t633 * t640;
	t673 = t578 * t628 - t591 * t653;
	t672 = t578 * t632 + t591 * t654;
	t660 = t625 * t631;
	t617 = -t627 * t634 + t630 * t660;
	t618 = t627 * t630 + t631 * t659;
	t601 = t617 * t629 + t618 * t633;
	t658 = t625 * t635;
	t649 = qJD(2) * t658;
	t606 = t618 * qJD(3) + t630 * t649;
	t607 = -t617 * qJD(3) + t634 * t649;
	t583 = t601 * qJD(5) - t606 * t633 + t607 * t629;
	t600 = t617 * t633 - t618 * t629;
	t671 = t583 * t628 - t600 * t653;
	t670 = t583 * t632 + t600 * t654;
	t663 = qJD(5) - qJD(3);
	t577 = t588 * qJD(5) + t593 * t629 + t594 * t633;
	t580 = t591 * qJD(5) + t595 * t629 + t596 * t633;
	t585 = t600 * qJD(5) + t606 * t629 + t607 * t633;
	t650 = t625 * t655;
	t642 = t629 * t634 - t630 * t633;
	t641 = t629 * t630 + t633 * t634;
	t637 = t663 * t642;
	t636 = t663 * t641;
	t613 = -t624 * t631 + t651;
	t612 = t638 * qJD(2);
	t610 = t614 * qJD(2);
	t608 = t641 * t658;
	t598 = t641 * t639;
	t597 = t641 * t613;
	t586 = (-t637 * t635 - t641 * t655) * t625;
	t582 = t641 * t612 + t637 * t639;
	t581 = -t641 * t610 - t637 * t613;
	t1 = [0, t582 * t632 + t611 * t628 + (t598 * t628 + t632 * t638) * qJD(6), t672, 0, -t672, -t580 * t628 + t612 * t632 + (-t592 * t632 + t628 * t639) * qJD(6); 0, t581 * t632 + t609 * t628 + (-t597 * t628 - t614 * t632) * qJD(6), t674, 0, -t674, -t577 * t628 - t610 * t632 + (-t589 * t632 - t613 * t628) * qJD(6); 0, -t628 * t649 + t586 * t632 + (-t608 * t628 - t632 * t660) * qJD(6), t670, 0, -t670, -t632 * t650 - t585 * t628 + (-t601 * t632 - t628 * t658) * qJD(6); 0, -t582 * t628 + t611 * t632 + (t598 * t632 - t628 * t638) * qJD(6), -t673, 0, t673, -t580 * t632 - t612 * t628 + (t592 * t628 + t632 * t639) * qJD(6); 0, -t581 * t628 + t609 * t632 + (-t597 * t632 + t614 * t628) * qJD(6), -t675, 0, t675, -t577 * t632 + t610 * t628 + (t589 * t628 - t613 * t632) * qJD(6); 0, -t632 * t649 - t586 * t628 + (-t608 * t632 + t628 * t660) * qJD(6), -t671, 0, t671, t628 * t650 - t585 * t632 + (t601 * t628 - t632 * t658) * qJD(6); 0, t642 * t612 - t636 * t639, -t580, 0, t580, 0; 0, -t642 * t610 + t636 * t613, -t577, 0, t577, 0; 0, (t636 * t635 - t642 * t655) * t625, -t585, 0, t585, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end