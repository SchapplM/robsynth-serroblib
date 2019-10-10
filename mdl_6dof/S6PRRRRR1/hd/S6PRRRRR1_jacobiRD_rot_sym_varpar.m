% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:13
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRRRR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:44
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:45
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
	t210 = sin(pkin(12));
	t212 = cos(pkin(12));
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
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:46
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (134->21), mult. (212->50), div. (0->0), fcn. (224->8), ass. (0->37)
	t274 = qJ(3) + qJ(4);
	t271 = sin(t274);
	t273 = qJD(3) + qJD(4);
	t293 = t271 * t273;
	t272 = cos(t274);
	t292 = t272 * t273;
	t276 = sin(pkin(6));
	t291 = t273 * t276;
	t280 = cos(qJ(2));
	t290 = t273 * t280;
	t278 = cos(pkin(6));
	t279 = sin(qJ(2));
	t289 = t278 * t279;
	t288 = t278 * t280;
	t287 = qJD(2) * t279;
	t286 = t279 * t291;
	t285 = t276 * qJD(2) * t280;
	t275 = sin(pkin(12));
	t277 = cos(pkin(12));
	t267 = -t275 * t279 + t277 * t288;
	t263 = t267 * qJD(2);
	t284 = t277 * t291 - t263;
	t269 = -t275 * t288 - t277 * t279;
	t265 = t269 * qJD(2);
	t283 = -t275 * t291 - t265;
	t268 = t275 * t280 + t277 * t289;
	t282 = t275 * t289 - t277 * t280;
	t281 = -t273 * t278 - t285;
	t266 = t282 * qJD(2);
	t264 = t268 * qJD(2);
	t262 = t271 * t286 + t281 * t272;
	t261 = t281 * t271 - t272 * t286;
	t260 = t283 * t272 - t282 * t293;
	t259 = t283 * t271 + t282 * t292;
	t258 = t268 * t293 + t284 * t272;
	t257 = -t268 * t292 + t284 * t271;
	t1 = [0, t266 * t272 - t269 * t293, t259, t259, 0, 0; 0, -t264 * t272 - t267 * t293, t257, t257, 0, 0; 0, (-t271 * t290 - t272 * t287) * t276, t261, t261, 0, 0; 0, -t266 * t271 - t269 * t292, t260, t260, 0, 0; 0, t264 * t271 - t267 * t292, t258, t258, 0, 0; 0, (t271 * t287 - t272 * t290) * t276, t262, t262, 0, 0; 0, t265, 0, 0, 0, 0; 0, t263, 0, 0, 0, 0; 0, t285, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:46
	% EndTime: 2019-10-09 23:13:46
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (292->21), mult. (284->50), div. (0->0), fcn. (300->8), ass. (0->37)
	t289 = qJ(3) + qJ(4) + qJ(5);
	t286 = sin(t289);
	t288 = qJD(3) + qJD(4) + qJD(5);
	t308 = t286 * t288;
	t287 = cos(t289);
	t307 = t287 * t288;
	t291 = sin(pkin(6));
	t306 = t288 * t291;
	t295 = cos(qJ(2));
	t305 = t288 * t295;
	t293 = cos(pkin(6));
	t294 = sin(qJ(2));
	t304 = t293 * t294;
	t303 = t293 * t295;
	t302 = qJD(2) * t294;
	t301 = t294 * t306;
	t300 = t291 * qJD(2) * t295;
	t290 = sin(pkin(12));
	t292 = cos(pkin(12));
	t282 = -t290 * t294 + t292 * t303;
	t278 = t282 * qJD(2);
	t299 = t292 * t306 - t278;
	t284 = -t290 * t303 - t292 * t294;
	t280 = t284 * qJD(2);
	t298 = -t290 * t306 - t280;
	t283 = t290 * t295 + t292 * t304;
	t297 = t290 * t304 - t292 * t295;
	t296 = -t288 * t293 - t300;
	t281 = t297 * qJD(2);
	t279 = t283 * qJD(2);
	t277 = t286 * t301 + t296 * t287;
	t276 = t296 * t286 - t287 * t301;
	t275 = t298 * t287 - t297 * t308;
	t274 = t298 * t286 + t297 * t307;
	t273 = t283 * t308 + t299 * t287;
	t272 = -t283 * t307 + t299 * t286;
	t1 = [0, t281 * t287 - t284 * t308, t274, t274, t274, 0; 0, -t279 * t287 - t282 * t308, t272, t272, t272, 0; 0, (-t286 * t305 - t287 * t302) * t291, t276, t276, t276, 0; 0, -t281 * t286 - t284 * t307, t275, t275, t275, 0; 0, t279 * t286 - t282 * t307, t273, t273, t273, 0; 0, (t286 * t302 - t287 * t305) * t291, t277, t277, t277, 0; 0, t280, 0, 0, 0, 0; 0, t278, 0, 0, 0, 0; 0, t300, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:48
	% EndTime: 2019-10-09 23:13:48
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (737->64), mult. (828->131), div. (0->0), fcn. (908->10), ass. (0->66)
	t461 = qJ(3) + qJ(4) + qJ(5);
	t458 = sin(t461);
	t459 = cos(t461);
	t467 = sin(qJ(2));
	t460 = qJD(3) + qJD(4) + qJD(5);
	t469 = cos(qJ(2));
	t496 = t460 * t469;
	t499 = (qJD(2) * t459 - qJD(6)) * t467 + t458 * t496;
	t498 = t458 * t460;
	t497 = t459 * t460;
	t462 = sin(pkin(12));
	t463 = sin(pkin(6));
	t495 = t462 * t463;
	t464 = cos(pkin(12));
	t494 = t463 * t464;
	t493 = t463 * t467;
	t492 = t463 * t469;
	t465 = cos(pkin(6));
	t491 = t465 * t467;
	t490 = t465 * t469;
	t489 = qJD(2) * t467;
	t488 = qJD(2) * t469;
	t487 = qJD(6) * t459;
	t466 = sin(qJ(6));
	t486 = qJD(6) * t466;
	t468 = cos(qJ(6));
	t485 = qJD(6) * t468;
	t483 = t462 * t491;
	t482 = t458 * t493;
	t481 = t459 * t493;
	t480 = t463 * t489;
	t473 = -t462 * t467 + t464 * t490;
	t448 = t473 * qJD(2);
	t478 = t460 * t494 - t448;
	t454 = t462 * t490 + t464 * t467;
	t450 = t454 * qJD(2);
	t477 = t460 * t495 - t450;
	t476 = -t473 * t487 + t448;
	t475 = t454 * t487 - t450;
	t474 = (qJD(2) - t487) * t469;
	t453 = t462 * t469 + t464 * t491;
	t472 = t460 * t465 + t463 * t488;
	t449 = t453 * qJD(2);
	t471 = qJD(6) * t453 - t449 * t459 - t473 * t498;
	t451 = -qJD(2) * t483 + t464 * t488;
	t455 = t464 * t469 - t483;
	t470 = qJD(6) * t455 - t451 * t459 + t454 * t498;
	t447 = t465 * t458 + t481;
	t446 = t465 * t459 - t482;
	t445 = t455 * t459 + t458 * t495;
	t444 = -t455 * t458 + t459 * t495;
	t443 = t453 * t459 - t458 * t494;
	t442 = -t453 * t458 - t459 * t494;
	t441 = t472 * t459 - t460 * t482;
	t440 = -t472 * t458 - t460 * t481;
	t439 = -t455 * t498 + t477 * t459;
	t438 = -t455 * t497 - t477 * t458;
	t437 = -t453 * t498 - t478 * t459;
	t436 = -t453 * t497 + t478 * t458;
	t435 = t440 * t468 - t446 * t486;
	t434 = -t440 * t466 - t446 * t485;
	t433 = t438 * t468 - t444 * t486;
	t432 = -t438 * t466 - t444 * t485;
	t431 = t436 * t468 - t442 * t486;
	t430 = -t436 * t466 - t442 * t485;
	t1 = [0, t475 * t466 + t470 * t468, t433, t433, t433, -t439 * t466 + t451 * t468 + (-t445 * t468 - t454 * t466) * qJD(6); 0, t476 * t466 + t471 * t468, t431, t431, t431, -t437 * t466 + t449 * t468 + (-t443 * t468 + t466 * t473) * qJD(6); 0, (t466 * t474 - t499 * t468) * t463, t435, t435, t435, t468 * t480 - t441 * t466 + (-t447 * t468 + t466 * t492) * qJD(6); 0, -t470 * t466 + t475 * t468, t432, t432, t432, -t439 * t468 - t451 * t466 + (t445 * t466 - t454 * t468) * qJD(6); 0, -t471 * t466 + t476 * t468, t430, t430, t430, -t437 * t468 - t449 * t466 + (t443 * t466 + t468 * t473) * qJD(6); 0, (t499 * t466 + t468 * t474) * t463, t434, t434, t434, -t466 * t480 - t441 * t468 + (t447 * t466 + t468 * t492) * qJD(6); 0, -t451 * t458 - t454 * t497, t439, t439, t439, 0; 0, -t449 * t458 + t473 * t497, t437, t437, t437, 0; 0, (-t458 * t489 + t459 * t496) * t463, t441, t441, t441, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end