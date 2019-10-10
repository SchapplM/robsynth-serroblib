% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:41
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPPRR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:14
	% EndTime: 2019-10-10 09:41:14
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (27->13), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t158 = sin(qJ(2));
	t159 = sin(qJ(1));
	t171 = t158 * t159;
	t161 = cos(qJ(1));
	t170 = t158 * t161;
	t160 = cos(qJ(2));
	t169 = t159 * t160;
	t168 = t160 * t161;
	t156 = sin(pkin(6));
	t167 = qJD(1) * t156;
	t166 = qJD(2) * t156;
	t157 = cos(pkin(6));
	t165 = -t157 * t168 + t171;
	t164 = t157 * t169 + t170;
	t163 = t157 * t170 + t169;
	t162 = t157 * t171 - t168;
	t155 = t162 * qJD(1) + t165 * qJD(2);
	t154 = t164 * qJD(1) + t163 * qJD(2);
	t153 = t163 * qJD(1) + t164 * qJD(2);
	t152 = t165 * qJD(1) + t162 * qJD(2);
	t1 = [t155, t152, 0, 0, 0, 0; -t153, -t154, 0, 0, 0, 0; 0, -t158 * t166, 0, 0, 0, 0; t154, t153, 0, 0, 0, 0; t152, t155, 0, 0, 0, 0; 0, -t160 * t166, 0, 0, 0, 0; -t159 * t167, 0, 0, 0, 0, 0; t161 * t167, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:14
	% EndTime: 2019-10-10 09:41:14
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (59->20), mult. (192->41), div. (0->0), fcn. (208->8), ass. (0->25)
	t224 = sin(pkin(11));
	t226 = cos(pkin(11));
	t227 = cos(pkin(6));
	t230 = cos(qJ(2));
	t237 = qJD(2) * t230;
	t228 = sin(qJ(2));
	t238 = qJD(2) * t228;
	t212 = (t224 * t238 - t226 * t237) * t227;
	t234 = t230 * t224 + t228 * t226;
	t215 = t234 * t227;
	t217 = -t224 * t237 - t226 * t238;
	t218 = t228 * t224 - t230 * t226;
	t229 = sin(qJ(1));
	t231 = cos(qJ(1));
	t240 = -t229 * t212 - t231 * t217 + (t215 * t231 - t218 * t229) * qJD(1);
	t225 = sin(pkin(6));
	t239 = qJD(1) * t225;
	t233 = qJD(2) * t234;
	t216 = t218 * qJD(2);
	t232 = t231 * t212 - t229 * t217 + (t215 * t229 + t218 * t231) * qJD(1);
	t214 = t218 * t227;
	t213 = t227 * t233;
	t211 = -t231 * t213 + t229 * t216 + (t214 * t229 - t231 * t234) * qJD(1);
	t210 = t229 * t213 + t231 * t216 + (t214 * t231 + t229 * t234) * qJD(1);
	t1 = [t232, t210, 0, 0, 0, 0; -t240, t211, 0, 0, 0, 0; 0, -t225 * t233, 0, 0, 0, 0; -t211, t240, 0, 0, 0, 0; t210, t232, 0, 0, 0, 0; 0, t225 * t216, 0, 0, 0, 0; -t229 * t239, 0, 0, 0, 0, 0; t231 * t239, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:14
	% EndTime: 2019-10-10 09:41:15
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (59->20), mult. (192->41), div. (0->0), fcn. (208->8), ass. (0->25)
	t282 = cos(pkin(6));
	t279 = sin(pkin(11));
	t281 = cos(pkin(11));
	t283 = sin(qJ(2));
	t285 = cos(qJ(2));
	t288 = t285 * t279 + t283 * t281;
	t270 = t288 * t282;
	t291 = qJD(2) * t285;
	t292 = qJD(2) * t283;
	t267 = (t279 * t292 - t281 * t291) * t282;
	t272 = -t279 * t291 - t281 * t292;
	t273 = t283 * t279 - t285 * t281;
	t284 = sin(qJ(1));
	t286 = cos(qJ(1));
	t295 = -t284 * t267 - t286 * t272 + (t270 * t286 - t273 * t284) * qJD(1);
	t280 = sin(pkin(6));
	t294 = qJD(1) * t280;
	t293 = qJD(2) * t280;
	t287 = -t286 * t267 + t284 * t272 + (-t270 * t284 - t273 * t286) * qJD(1);
	t271 = t273 * qJD(2);
	t269 = t273 * t282;
	t268 = qJD(2) * t270;
	t266 = -t286 * t268 + t284 * t271 + (t269 * t284 - t286 * t288) * qJD(1);
	t265 = -t284 * t268 - t286 * t271 + (-t269 * t286 - t284 * t288) * qJD(1);
	t1 = [-t284 * t294, 0, 0, 0, 0, 0; t286 * t294, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t287, t265, 0, 0, 0, 0; t295, -t266, 0, 0, 0, 0; 0, t288 * t293, 0, 0, 0, 0; t266, -t295, 0, 0, 0, 0; t265, t287, 0, 0, 0, 0; 0, -t273 * t293, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:15
	% EndTime: 2019-10-10 09:41:16
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (191->43), mult. (594->87), div. (0->0), fcn. (666->10), ass. (0->44)
	t385 = cos(pkin(6));
	t382 = sin(pkin(11));
	t384 = cos(pkin(11));
	t387 = sin(qJ(2));
	t390 = cos(qJ(2));
	t395 = t390 * t382 + t387 * t384;
	t373 = t395 * t385;
	t376 = t387 * t382 - t390 * t384;
	t374 = t376 * qJD(2);
	t402 = qJD(2) * t390;
	t403 = qJD(2) * t387;
	t368 = (t382 * t403 - t384 * t402) * t385;
	t375 = -t382 * t402 - t384 * t403;
	t388 = sin(qJ(1));
	t391 = cos(qJ(1));
	t396 = t388 * t373 + t391 * t376;
	t408 = t396 * qJD(1) + t391 * t368 - t388 * t375;
	t369 = qJD(2) * t373;
	t372 = t376 * t385;
	t404 = qJD(1) * t388;
	t360 = t388 * t374 + t372 * t404 + (-qJD(1) * t395 - t369) * t391;
	t362 = -t391 * t372 - t388 * t395;
	t386 = sin(qJ(5));
	t389 = cos(qJ(5));
	t383 = sin(pkin(6));
	t399 = t383 * t404;
	t405 = t383 * t391;
	t407 = (t362 * t386 + t389 * t405) * qJD(5) - t360 * t389 - t386 * t399;
	t406 = t383 * t388;
	t401 = qJD(5) * t386;
	t400 = qJD(5) * t389;
	t398 = qJD(1) * t405;
	t361 = t391 * t373 - t388 * t376;
	t371 = t395 * t383;
	t357 = -t361 * qJD(1) + t388 * t368 + t391 * t375;
	t392 = -t389 * t399 + t360 * t386 + (t362 * t389 - t386 * t405) * qJD(5);
	t370 = t376 * t383;
	t367 = t383 * t374;
	t366 = qJD(2) * t371;
	t364 = -t388 * t372 + t391 * t395;
	t358 = t362 * qJD(1) - t388 * t369 - t391 * t374;
	t356 = t389 * t398 + t358 * t386 + (t364 * t389 - t386 * t406) * qJD(5);
	t355 = -t386 * t398 + t358 * t389 + (-t364 * t386 - t389 * t406) * qJD(5);
	t1 = [t392, t357 * t386 - t396 * t400, 0, 0, t355, 0; t356, t361 * t400 - t386 * t408, 0, 0, t407, 0; 0, -t367 * t386 + t371 * t400, 0, 0, t366 * t389 + (-t370 * t386 - t385 * t389) * qJD(5), 0; -t407, t357 * t389 + t396 * t401, 0, 0, -t356, 0; t355, -t361 * t401 - t389 * t408, 0, 0, t392, 0; 0, -t367 * t389 - t371 * t401, 0, 0, -t366 * t386 + (-t370 * t389 + t385 * t386) * qJD(5), 0; t408, -t358, 0, 0, 0, 0; t357, t360, 0, 0, 0, 0; 0, -t366, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:18
	% EndTime: 2019-10-10 09:41:19
	% DurationCPUTime: 0.71s
	% Computational Cost: add. (555->81), mult. (1658->154), div. (0->0), fcn. (1910->12), ass. (0->71)
	t596 = cos(pkin(6));
	t593 = sin(pkin(11));
	t595 = cos(pkin(11));
	t599 = sin(qJ(2));
	t603 = cos(qJ(2));
	t613 = t603 * t593 + t599 * t595;
	t583 = t613 * t596;
	t579 = qJD(2) * t583;
	t600 = sin(qJ(1));
	t604 = cos(qJ(1));
	t585 = t599 * t593 - t603 * t595;
	t612 = t585 * t596;
	t609 = t600 * t612;
	t611 = t585 * qJD(2);
	t556 = qJD(1) * t609 + (-qJD(1) * t613 - t579) * t604 + t600 * t611;
	t566 = -t600 * t613 - t604 * t612;
	t598 = sin(qJ(5));
	t602 = cos(qJ(5));
	t594 = sin(pkin(6));
	t633 = t594 * t604;
	t562 = -t566 * t602 + t598 * t633;
	t629 = qJD(1) * t600;
	t620 = t594 * t629;
	t547 = t562 * qJD(5) - t556 * t598 + t602 * t620;
	t627 = qJD(2) * t603;
	t628 = qJD(2) * t599;
	t578 = (t593 * t628 - t595 * t627) * t596;
	t584 = -t593 * t627 - t595 * t628;
	t631 = t600 * t584;
	t557 = -t583 * t629 + t631 + (-qJD(1) * t585 - t578) * t604;
	t621 = t602 * t633;
	t563 = t566 * t598 + t621;
	t597 = sin(qJ(6));
	t601 = cos(qJ(6));
	t615 = t604 * t583 - t600 * t585;
	t642 = -t547 * t601 + (-t563 * t597 - t615 * t601) * qJD(6) - t557 * t597;
	t641 = (t563 * t601 - t615 * t597) * qJD(6) - t547 * t597 + t557 * t601;
	t634 = t594 * t600;
	t626 = qJD(5) * t598;
	t625 = qJD(5) * t602;
	t624 = qJD(6) * t597;
	t623 = qJD(6) * t598;
	t622 = qJD(6) * t601;
	t619 = qJD(1) * t633;
	t553 = t566 * qJD(1) - t600 * t579 - t604 * t611;
	t614 = -t600 * t583 - t604 * t585;
	t618 = -t614 * t623 - t553;
	t617 = -t615 * t623 + t556;
	t582 = t613 * t594;
	t576 = qJD(2) * t582;
	t616 = -t582 * t623 - t576;
	t581 = t585 * t594;
	t571 = t581 * t602 - t596 * t598;
	t572 = t581 * t598 + t596 * t602;
	t569 = -t604 * t613 + t609;
	t561 = -t569 * t598 + t602 * t634;
	t560 = -t569 * t602 - t598 * t634;
	t605 = -t615 * qJD(1) + t600 * t578 + t604 * t584;
	t608 = qJD(6) * t569 + t598 * t605 + t614 * t625;
	t555 = t614 * qJD(1) - t604 * t578 + t631;
	t607 = qJD(6) * t566 + t555 * t598 + t615 * t625;
	t577 = t594 * t611;
	t606 = -qJD(6) * t581 - t577 * t598 + t582 * t625;
	t546 = -t556 * t602 + qJD(5) * t621 + (qJD(5) * t566 - t620) * t598;
	t559 = -t572 * qJD(5) + t576 * t602;
	t558 = t571 * qJD(5) + t576 * t598;
	t550 = t560 * qJD(5) + t553 * t598 + t602 * t619;
	t549 = t561 * qJD(5) - t553 * t602 + t598 * t619;
	t545 = t550 * t601 + t605 * t597 + (-t561 * t597 + t601 * t614) * qJD(6);
	t544 = -t550 * t597 + t605 * t601 + (-t561 * t601 - t597 * t614) * qJD(6);
	t1 = [t642, t618 * t597 + t608 * t601, 0, 0, -t549 * t601 - t560 * t624, t544; t545, t617 * t597 + t607 * t601, 0, 0, t546 * t601 - t562 * t624, t641; 0, t616 * t597 + t606 * t601, 0, 0, t559 * t601 - t571 * t624, -t558 * t597 - t577 * t601 + (-t572 * t601 - t582 * t597) * qJD(6); -t641, -t608 * t597 + t618 * t601, 0, 0, t549 * t597 - t560 * t622, -t545; t544, -t607 * t597 + t617 * t601, 0, 0, -t546 * t597 - t562 * t622, t642; 0, -t606 * t597 + t616 * t601, 0, 0, -t559 * t597 - t571 * t622, -t558 * t601 + t577 * t597 + (t572 * t597 - t582 * t601) * qJD(6); t546, -t602 * t605 + t614 * t626, 0, 0, t550, 0; t549, -t555 * t602 + t615 * t626, 0, 0, t547, 0; 0, t577 * t602 + t582 * t626, 0, 0, t558, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end