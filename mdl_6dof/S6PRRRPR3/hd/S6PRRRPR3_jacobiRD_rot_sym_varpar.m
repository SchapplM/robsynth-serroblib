% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:50
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRRPR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:24
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
	% StartTime: 2019-10-09 22:50:25
	% EndTime: 2019-10-09 22:50:25
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
	% StartTime: 2019-10-09 22:50:25
	% EndTime: 2019-10-09 22:50:25
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
	t275 = sin(pkin(11));
	t277 = cos(pkin(11));
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
	% StartTime: 2019-10-09 22:50:26
	% EndTime: 2019-10-09 22:50:26
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (134->21), mult. (212->50), div. (0->0), fcn. (224->8), ass. (0->37)
	t341 = qJ(3) + qJ(4);
	t338 = sin(t341);
	t340 = qJD(3) + qJD(4);
	t360 = t338 * t340;
	t339 = cos(t341);
	t359 = t339 * t340;
	t343 = sin(pkin(6));
	t358 = t340 * t343;
	t347 = cos(qJ(2));
	t357 = t340 * t347;
	t345 = cos(pkin(6));
	t346 = sin(qJ(2));
	t356 = t345 * t346;
	t355 = t345 * t347;
	t354 = qJD(2) * t346;
	t353 = t346 * t358;
	t352 = t343 * qJD(2) * t347;
	t342 = sin(pkin(11));
	t344 = cos(pkin(11));
	t334 = -t342 * t346 + t344 * t355;
	t330 = t334 * qJD(2);
	t351 = -t344 * t358 + t330;
	t336 = -t342 * t355 - t344 * t346;
	t332 = t336 * qJD(2);
	t350 = t342 * t358 + t332;
	t335 = t342 * t347 + t344 * t356;
	t349 = t342 * t356 - t344 * t347;
	t348 = t340 * t345 + t352;
	t333 = t349 * qJD(2);
	t331 = t335 * qJD(2);
	t329 = -t338 * t353 + t348 * t339;
	t328 = t348 * t338 + t339 * t353;
	t327 = t350 * t339 + t349 * t360;
	t326 = t350 * t338 - t349 * t359;
	t325 = -t335 * t360 + t351 * t339;
	t324 = t335 * t359 + t351 * t338;
	t1 = [0, t332, 0, 0, 0, 0; 0, t330, 0, 0, 0, 0; 0, t352, 0, 0, 0, 0; 0, -t333 * t339 + t336 * t360, t326, t326, 0, 0; 0, t331 * t339 + t334 * t360, t324, t324, 0, 0; 0, (t338 * t357 + t339 * t354) * t343, t328, t328, 0, 0; 0, t333 * t338 + t336 * t359, t327, t327, 0, 0; 0, -t331 * t338 + t334 * t359, t325, t325, 0, 0; 0, (-t338 * t354 + t339 * t357) * t343, t329, t329, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:50:27
	% EndTime: 2019-10-09 22:50:28
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (388->71), mult. (672->132), div. (0->0), fcn. (736->10), ass. (0->67)
	t467 = qJ(3) + qJ(4);
	t464 = sin(t467);
	t466 = qJD(3) + qJD(4);
	t504 = t464 * t466;
	t465 = cos(t467);
	t503 = t465 * t466;
	t475 = cos(qJ(2));
	t502 = t466 * t475;
	t468 = sin(pkin(11));
	t469 = sin(pkin(6));
	t501 = t468 * t469;
	t470 = cos(pkin(11));
	t500 = t469 * t470;
	t473 = sin(qJ(2));
	t499 = t469 * t473;
	t471 = cos(pkin(6));
	t498 = t471 * t473;
	t497 = t471 * t475;
	t472 = sin(qJ(6));
	t496 = t472 * t475;
	t474 = cos(qJ(6));
	t495 = t474 * t475;
	t494 = qJD(2) * t473;
	t493 = qJD(2) * t475;
	t492 = qJD(6) * t464;
	t491 = qJD(6) * t472;
	t490 = qJD(6) * t474;
	t489 = t468 * t498;
	t488 = t464 * t499;
	t487 = t465 * t499;
	t486 = t464 * t500;
	t485 = t469 * t494;
	t484 = qJD(2) + t492;
	t459 = t468 * t497 + t470 * t473;
	t455 = t459 * qJD(2);
	t483 = t466 * t501 - t455;
	t480 = -t468 * t473 + t470 * t497;
	t453 = t480 * qJD(2);
	t482 = -t480 * t492 - t453;
	t481 = t459 * t492 + t455;
	t458 = t468 * t475 + t470 * t498;
	t479 = t466 * t471 + t469 * t493;
	t454 = t458 * qJD(2);
	t478 = -qJD(6) * t458 - t454 * t464 + t480 * t503;
	t456 = -qJD(2) * t489 + t470 * t493;
	t460 = t470 * t475 - t489;
	t477 = -qJD(6) * t460 - t456 * t464 - t459 * t503;
	t476 = t465 * t502 + (-qJD(2) * t464 - qJD(6)) * t473;
	t452 = t471 * t464 + t487;
	t451 = -t471 * t465 + t488;
	t450 = t460 * t465 + t464 * t501;
	t449 = t460 * t464 - t465 * t501;
	t448 = t458 * t465 - t486;
	t447 = t458 * t464 + t465 * t500;
	t446 = t479 * t465 - t466 * t488;
	t445 = t479 * t464 + t466 * t487;
	t444 = -t460 * t504 + t483 * t465;
	t443 = t460 * t503 + t483 * t464;
	t442 = -t458 * t504 + (-t466 * t500 + t453) * t465;
	t441 = t453 * t464 + t458 * t503 - t466 * t486;
	t440 = t446 * t474 - t452 * t491;
	t439 = t446 * t472 + t452 * t490;
	t438 = t444 * t474 - t450 * t491;
	t437 = t444 * t472 + t450 * t490;
	t436 = t442 * t474 - t448 * t491;
	t435 = t442 * t472 + t448 * t490;
	t1 = [0, t477 * t472 - t481 * t474, t437, t437, 0, t443 * t474 - t456 * t472 + (-t449 * t472 - t459 * t474) * qJD(6); 0, t478 * t472 - t482 * t474, t435, t435, 0, t441 * t474 - t454 * t472 + (-t447 * t472 + t474 * t480) * qJD(6); 0, (t476 * t472 + t484 * t495) * t469, t439, t439, 0, -t472 * t485 + t445 * t474 + (-t451 * t472 + t469 * t495) * qJD(6); 0, t481 * t472 + t477 * t474, t438, t438, 0, -t443 * t472 - t456 * t474 + (-t449 * t474 + t459 * t472) * qJD(6); 0, t482 * t472 + t478 * t474, t436, t436, 0, -t441 * t472 - t454 * t474 + (-t447 * t474 - t472 * t480) * qJD(6); 0, (t476 * t474 - t484 * t496) * t469, t440, t440, 0, -t474 * t485 - t445 * t472 + (-t451 * t474 - t469 * t496) * qJD(6); 0, -t456 * t465 + t459 * t504, -t443, -t443, 0, 0; 0, -t454 * t465 - t480 * t504, -t441, -t441, 0, 0; 0, (-t464 * t502 - t465 * t494) * t469, -t445, -t445, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end