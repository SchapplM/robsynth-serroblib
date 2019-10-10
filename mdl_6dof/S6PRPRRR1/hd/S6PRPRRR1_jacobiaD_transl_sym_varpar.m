% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRRR1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:53
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(11));
	t50 = sin(pkin(11));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (22->14), mult. (83->33), div. (0->0), fcn. (72->8), ass. (0->16)
	t89 = sin(pkin(12));
	t92 = cos(pkin(12));
	t95 = sin(qJ(2));
	t96 = cos(qJ(2));
	t98 = t89 * t96 + t92 * t95;
	t88 = t98 * qJD(2);
	t94 = cos(pkin(6));
	t100 = t94 * t95;
	t99 = pkin(2) * qJD(2);
	t97 = t89 * t95 - t92 * t96;
	t87 = t97 * qJD(2);
	t93 = cos(pkin(11));
	t90 = sin(pkin(11));
	t86 = t94 * t88;
	t85 = t94 * t87;
	t1 = [0, (t90 * t86 + t93 * t87) * r_i_i_C(1) + (-t90 * t85 + t93 * t88) * r_i_i_C(2) + (t90 * t100 - t93 * t96) * t99, 0, 0, 0, 0; 0, (-t93 * t86 + t90 * t87) * r_i_i_C(1) + (t93 * t85 + t90 * t88) * r_i_i_C(2) + (-t93 * t100 - t90 * t96) * t99, 0, 0, 0, 0; 0, (-t95 * pkin(2) - t98 * r_i_i_C(1) + t97 * r_i_i_C(2)) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:49
	% EndTime: 2019-10-09 21:53:49
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (115->39), mult. (386->84), div. (0->0), fcn. (382->10), ass. (0->37)
	t241 = sin(pkin(12));
	t244 = cos(pkin(12));
	t250 = cos(qJ(2));
	t261 = qJD(2) * t250;
	t248 = sin(qJ(2));
	t262 = qJD(2) * t248;
	t268 = t241 * t262 - t244 * t261;
	t267 = -pkin(8) - r_i_i_C(3);
	t266 = pkin(2) * qJD(2);
	t243 = sin(pkin(6));
	t247 = sin(qJ(4));
	t265 = t243 * t247;
	t249 = cos(qJ(4));
	t264 = t243 * t249;
	t246 = cos(pkin(6));
	t263 = t246 * t248;
	t258 = r_i_i_C(1) * t247 + r_i_i_C(2) * t249;
	t227 = t268 * t246;
	t234 = -t241 * t261 - t244 * t262;
	t242 = sin(pkin(11));
	t245 = cos(pkin(11));
	t257 = t227 * t245 - t234 * t242;
	t256 = t227 * t242 + t234 * t245;
	t255 = t241 * t250 + t248 * t244;
	t254 = t248 * t241 - t244 * t250;
	t253 = r_i_i_C(1) * t249 - r_i_i_C(2) * t247 + pkin(3);
	t252 = qJD(4) * t258;
	t251 = qJD(2) * t255;
	t233 = t254 * qJD(2);
	t232 = t255 * t246;
	t231 = t254 * t246;
	t230 = t255 * t243;
	t228 = t246 * t251;
	t225 = t268 * t243;
	t224 = -t232 * t242 - t245 * t254;
	t222 = t232 * t245 - t242 * t254;
	t1 = [0, -t267 * t256 - (t231 * t242 - t245 * t255) * t252 + (t242 * t263 - t245 * t250) * t266 + t253 * (t228 * t242 + t233 * t245), 0, -t258 * t256 + ((-t224 * t249 - t242 * t265) * r_i_i_C(1) + (t224 * t247 - t242 * t264) * r_i_i_C(2)) * qJD(4), 0, 0; 0, t267 * t257 - (-t231 * t245 - t242 * t255) * t252 + (-t242 * t250 - t245 * t263) * t266 + t253 * (-t228 * t245 + t233 * t242), 0, t258 * t257 + ((-t222 * t249 + t245 * t265) * r_i_i_C(1) + (t222 * t247 + t245 * t264) * r_i_i_C(2)) * qJD(4), 0, 0; 0, t267 * t225 + (-pkin(2) * t262 - t251 * t253 + t252 * t254) * t243, 0, t258 * t225 + ((-t230 * t249 - t246 * t247) * r_i_i_C(1) + (t230 * t247 - t246 * t249) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:49
	% EndTime: 2019-10-09 21:53:49
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (275->52), mult. (634->101), div. (0->0), fcn. (640->12), ass. (0->48)
	t283 = sin(pkin(12));
	t286 = cos(pkin(12));
	t292 = cos(qJ(2));
	t305 = qJD(2) * t292;
	t290 = sin(qJ(2));
	t306 = qJD(2) * t290;
	t317 = t283 * t306 - t286 * t305;
	t316 = -r_i_i_C(3) - pkin(9) - pkin(8);
	t315 = pkin(2) * qJD(2);
	t282 = qJ(4) + qJ(5);
	t279 = sin(t282);
	t281 = qJD(4) + qJD(5);
	t314 = t279 * t281;
	t280 = cos(t282);
	t313 = t280 * t281;
	t285 = sin(pkin(6));
	t312 = t281 * t285;
	t289 = sin(qJ(4));
	t311 = t285 * t289;
	t288 = cos(pkin(6));
	t310 = t288 * t290;
	t298 = t292 * t283 + t290 * t286;
	t269 = t298 * t288;
	t284 = sin(pkin(11));
	t287 = cos(pkin(11));
	t297 = t290 * t283 - t292 * t286;
	t259 = t287 * t269 - t284 * t297;
	t264 = t317 * t288;
	t271 = -t283 * t305 - t286 * t306;
	t299 = t287 * t264 - t284 * t271;
	t301 = t287 * t312 + t299;
	t309 = (-t259 * t313 + t301 * t279) * r_i_i_C(1) + (t259 * t314 + t301 * t280) * r_i_i_C(2);
	t261 = -t284 * t269 - t287 * t297;
	t256 = t284 * t264 + t287 * t271;
	t300 = -t284 * t312 - t256;
	t308 = (-t261 * t313 + t300 * t279) * r_i_i_C(1) + (t261 * t314 + t300 * t280) * r_i_i_C(2);
	t267 = t298 * t285;
	t262 = t317 * t285;
	t302 = -t281 * t288 + t262;
	t307 = (-t267 * t313 + t302 * t279) * r_i_i_C(1) + (t267 * t314 + t302 * t280) * r_i_i_C(2);
	t291 = cos(qJ(4));
	t296 = t291 * pkin(4) + r_i_i_C(1) * t280 - r_i_i_C(2) * t279 + pkin(3);
	t295 = qJD(2) * t298;
	t294 = -pkin(4) * qJD(4) * t289 + (-r_i_i_C(1) * t279 - r_i_i_C(2) * t280) * t281;
	t270 = t297 * qJD(2);
	t268 = t297 * t288;
	t265 = t288 * t295;
	t1 = [0, -t316 * t256 + (t284 * t310 - t287 * t292) * t315 + t294 * (t284 * t268 - t287 * t298) + t296 * (t284 * t265 + t287 * t270), 0, (-t256 * t289 + (-t261 * t291 - t284 * t311) * qJD(4)) * pkin(4) + t308, t308, 0; 0, t316 * t299 + (-t284 * t292 - t287 * t310) * t315 + t294 * (-t287 * t268 - t284 * t298) + t296 * (-t287 * t265 + t284 * t270), 0, (t299 * t289 + (-t259 * t291 + t287 * t311) * qJD(4)) * pkin(4) + t309, t309, 0; 0, t316 * t262 + (-pkin(2) * t306 - t294 * t297 - t295 * t296) * t285, 0, (t262 * t289 + (-t267 * t291 - t288 * t289) * qJD(4)) * pkin(4) + t307, t307, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:51
	% EndTime: 2019-10-09 21:53:52
	% DurationCPUTime: 0.92s
	% Computational Cost: add. (873->100), mult. (1890->184), div. (0->0), fcn. (2030->14), ass. (0->67)
	t490 = pkin(10) + r_i_i_C(3);
	t449 = cos(qJ(6));
	t475 = qJD(6) * t449;
	t446 = sin(qJ(6));
	t476 = qJD(6) * t446;
	t495 = -r_i_i_C(1) * t476 - t475 * r_i_i_C(2);
	t440 = qJ(4) + qJ(5);
	t437 = sin(t440);
	t438 = cos(t440);
	t439 = qJD(4) + qJD(5);
	t447 = sin(qJ(4));
	t462 = t449 * r_i_i_C(1) - t446 * r_i_i_C(2) + pkin(5);
	t494 = (t462 * t437 - t490 * t438) * t439 + (t446 * r_i_i_C(1) + t449 * r_i_i_C(2)) * t438 * qJD(6) + qJD(4) * t447 * pkin(4);
	t445 = cos(pkin(6));
	t441 = sin(pkin(12));
	t448 = sin(qJ(2));
	t486 = cos(pkin(12));
	t489 = cos(qJ(2));
	t460 = t489 * t441 + t448 * t486;
	t424 = t460 * t445;
	t467 = t489 * t486;
	t478 = qJD(2) * t448;
	t493 = -qJD(2) * t467 + t441 * t478;
	t459 = -t448 * t441 + t467;
	t450 = cos(qJ(4));
	t456 = t450 * pkin(4) + t490 * t437 + t462 * t438 + pkin(3);
	t487 = pkin(2) * qJD(2);
	t485 = t437 * t439;
	t484 = t438 * t439;
	t442 = sin(pkin(11));
	t443 = sin(pkin(6));
	t483 = t442 * t443;
	t444 = cos(pkin(11));
	t482 = t443 * t444;
	t481 = t443 * t447;
	t480 = t445 * t448;
	t472 = t438 * t482;
	t419 = t493 * t443;
	t468 = t439 * t445 - t419;
	t421 = t493 * t445;
	t426 = t460 * qJD(2);
	t405 = t442 * t421 - t444 * t426;
	t465 = t439 * t483 + t405;
	t403 = t444 * t421 + t442 * t426;
	t464 = t444 * t424 + t442 * t459;
	t463 = -t442 * t424 + t444 * t459;
	t423 = t460 * t443;
	t458 = t459 * t445;
	t457 = qJD(2) * t424;
	t387 = -t403 * t438 - t439 * t472 - t464 * t485;
	t455 = t495 * (-t437 * t464 - t472) + t490 * t387 + t462 * (-t464 * t484 + (t439 * t482 + t403) * t437);
	t389 = t465 * t438 - t463 * t485;
	t454 = t495 * (-t437 * t463 + t438 * t483) + t490 * t389 + t462 * (-t465 * t437 - t463 * t484);
	t394 = -t423 * t485 + t468 * t438;
	t453 = t495 * (-t423 * t437 + t445 * t438) + t490 * t394 + t462 * (-t423 * t484 - t468 * t437);
	t451 = -pkin(9) - pkin(8);
	t425 = t459 * qJD(2);
	t422 = t459 * t443;
	t420 = qJD(2) * t423;
	t414 = t423 * t438 + t445 * t437;
	t411 = -t442 * t458 - t444 * t460;
	t408 = -t442 * t460 + t444 * t458;
	t404 = -t444 * t425 + t442 * t457;
	t401 = -t442 * t425 - t444 * t457;
	t398 = t437 * t483 + t438 * t463;
	t396 = -t437 * t482 + t438 * t464;
	t1 = [0, (t405 * t446 + t463 * t475) * r_i_i_C(1) + (t405 * t449 - t463 * t476) * r_i_i_C(2) - t405 * t451 + (t442 * t480 - t489 * t444) * t487 + t456 * t404 - t494 * t411, 0, (-t405 * t447 + (-t442 * t481 - t450 * t463) * qJD(4)) * pkin(4) + t454, t454, (-t389 * t446 - t404 * t449) * r_i_i_C(1) + (-t389 * t449 + t404 * t446) * r_i_i_C(2) + ((-t398 * t449 + t411 * t446) * r_i_i_C(1) + (t398 * t446 + t411 * t449) * r_i_i_C(2)) * qJD(6); 0, (-t403 * t446 + t464 * t475) * r_i_i_C(1) + (-t403 * t449 - t464 * t476) * r_i_i_C(2) + t403 * t451 + (-t489 * t442 - t444 * t480) * t487 + t456 * t401 - t494 * t408, 0, (t403 * t447 + (t444 * t481 - t450 * t464) * qJD(4)) * pkin(4) + t455, t455, (-t387 * t446 - t401 * t449) * r_i_i_C(1) + (-t387 * t449 + t401 * t446) * r_i_i_C(2) + ((-t396 * t449 + t408 * t446) * r_i_i_C(1) + (t396 * t446 + t408 * t449) * r_i_i_C(2)) * qJD(6); 0, (-t419 * t446 + t423 * t475) * r_i_i_C(1) + (-t419 * t449 - t423 * t476) * r_i_i_C(2) + t419 * t451 - t443 * pkin(2) * t478 - t456 * t420 - t494 * t422, 0, (t419 * t447 + (-t423 * t450 - t445 * t447) * qJD(4)) * pkin(4) + t453, t453, (-t394 * t446 + t420 * t449) * r_i_i_C(1) + (-t394 * t449 - t420 * t446) * r_i_i_C(2) + ((-t414 * t449 + t422 * t446) * r_i_i_C(1) + (t414 * t446 + t422 * t449) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end