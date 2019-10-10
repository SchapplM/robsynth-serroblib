% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR12
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:07
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR12_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR12_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR12_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (35->18), mult. (110->35), div. (0->0), fcn. (94->6), ass. (0->20)
	t136 = sin(pkin(6));
	t151 = t136 * (pkin(8) + r_i_i_C(3));
	t138 = sin(qJ(2));
	t139 = sin(qJ(1));
	t149 = t138 * t139;
	t141 = cos(qJ(1));
	t148 = t138 * t141;
	t140 = cos(qJ(2));
	t147 = t139 * t140;
	t146 = t140 * t141;
	t137 = cos(pkin(6));
	t145 = -t137 * t146 + t149;
	t144 = t137 * t147 + t148;
	t143 = t137 * t148 + t147;
	t142 = t137 * t149 - t146;
	t135 = t142 * qJD(1) + t145 * qJD(2);
	t134 = t144 * qJD(1) + t143 * qJD(2);
	t133 = t143 * qJD(1) + t144 * qJD(2);
	t132 = t145 * qJD(1) + t142 * qJD(2);
	t1 = [t135 * r_i_i_C(1) + t134 * r_i_i_C(2) + (-pkin(1) * t141 - t139 * t151) * qJD(1), t132 * r_i_i_C(1) + t133 * r_i_i_C(2), 0, 0, 0, 0; -t133 * r_i_i_C(1) + t132 * r_i_i_C(2) + (-pkin(1) * t139 + t141 * t151) * qJD(1), -t134 * r_i_i_C(1) + t135 * r_i_i_C(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t138 - r_i_i_C(2) * t140) * t136 * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:03
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (84->29), mult. (244->45), div. (0->0), fcn. (222->6), ass. (0->24)
	t175 = sin(pkin(6));
	t194 = t175 * (pkin(8) + r_i_i_C(1));
	t193 = pkin(2) - r_i_i_C(2);
	t191 = r_i_i_C(3) + qJ(3);
	t177 = sin(qJ(2));
	t178 = sin(qJ(1));
	t190 = t178 * t177;
	t179 = cos(qJ(2));
	t189 = t178 * t179;
	t180 = cos(qJ(1));
	t188 = t180 * t177;
	t187 = t180 * t179;
	t186 = qJD(2) * t177;
	t176 = cos(pkin(6));
	t185 = t176 * t190;
	t184 = t176 * t187;
	t183 = qJD(2) * t176 + qJD(1);
	t182 = t176 * t189 + t188;
	t181 = t176 * t188 + t189;
	t170 = -qJD(1) * t185 - t178 * t186 + t183 * t187;
	t169 = t182 * qJD(1) + t181 * qJD(2);
	t168 = t181 * qJD(1) + t182 * qJD(2);
	t167 = -qJD(1) * t184 - qJD(2) * t187 + t183 * t190;
	t1 = [-(-t184 + t190) * qJD(3) - t193 * t170 - t191 * t169 + (-t180 * pkin(1) - t178 * t194) * qJD(1), -(t185 - t187) * qJD(3) - t191 * t168 + t193 * t167, -t167, 0, 0, 0; t182 * qJD(3) - t193 * t168 - t191 * t167 + (-t178 * pkin(1) + t180 * t194) * qJD(1), t181 * qJD(3) - t193 * t169 + t191 * t170, t169, 0, 0, 0; 0, (t177 * qJD(3) + (-t193 * t177 + t191 * t179) * qJD(2)) * t175, t175 * t186, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:03
	% EndTime: 2019-10-10 11:07:03
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (179->58), mult. (534->98), div. (0->0), fcn. (502->8), ass. (0->43)
	t275 = pkin(3) + pkin(8);
	t243 = cos(pkin(6));
	t245 = sin(qJ(2));
	t249 = cos(qJ(1));
	t267 = t249 * t245;
	t246 = sin(qJ(1));
	t248 = cos(qJ(2));
	t268 = t246 * t248;
	t235 = t243 * t268 + t267;
	t251 = t243 * t267 + t268;
	t231 = t235 * qJD(1) + t251 * qJD(2);
	t244 = sin(qJ(4));
	t274 = t231 * t244;
	t247 = cos(qJ(4));
	t273 = t231 * t247;
	t242 = sin(pkin(6));
	t272 = t242 * t246;
	t271 = t242 * t248;
	t270 = t242 * t249;
	t269 = t246 * t245;
	t266 = t249 * t248;
	t265 = qJD(1) * t246;
	t264 = qJD(1) * t249;
	t263 = qJD(2) * t245;
	t259 = t243 * t266;
	t233 = -t259 + t269;
	t262 = qJD(4) * t233;
	t261 = -r_i_i_C(3) - pkin(9) - pkin(2);
	t260 = t243 * t269;
	t258 = t242 * t265;
	t257 = t242 * t264;
	t256 = t242 * t263;
	t255 = qJD(2) * t243 + qJD(1);
	t254 = r_i_i_C(1) * t247 - r_i_i_C(2) * t244;
	t253 = -t244 * r_i_i_C(1) - t247 * r_i_i_C(2);
	t252 = qJ(3) - t253;
	t250 = t254 * qJD(4) + qJD(3);
	t232 = -qJD(1) * t260 - t246 * t263 + t255 * t266;
	t230 = t251 * qJD(1) + t235 * qJD(2);
	t229 = -qJD(1) * t259 - qJD(2) * t266 + t255 * t269;
	t228 = t247 * t257 - t229 * t244 + (t235 * t247 - t244 * t272) * qJD(4);
	t227 = -t244 * t257 - t229 * t247 + (-t235 * t244 - t247 * t272) * qJD(4);
	t1 = [(-t247 * t262 - t274) * r_i_i_C(1) + (t244 * t262 - t273) * r_i_i_C(2) - t231 * qJ(3) - t233 * qJD(3) - pkin(1) * t264 + t261 * t232 + (t253 * t249 * qJD(4) + (-t254 - t275) * t265) * t242, t250 * (-t260 + t266) - t252 * t230 - t261 * t229, -t229, t227 * r_i_i_C(1) - t228 * r_i_i_C(2), 0, 0; t228 * r_i_i_C(1) + t227 * r_i_i_C(2) - t229 * qJ(3) + t235 * qJD(3) + t261 * t230 + (-pkin(1) * t246 + t275 * t270) * qJD(1), t261 * t231 + t252 * t232 + t250 * t251, t231, (-t244 * t258 + t273) * r_i_i_C(1) + (-t247 * t258 - t274) * r_i_i_C(2) + ((-t233 * t244 + t247 * t270) * r_i_i_C(1) + (-t233 * t247 - t244 * t270) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (t250 * t245 + (t261 * t245 + t252 * t248) * qJD(2)) * t242, t256, t254 * t256 + ((-t243 * t247 + t244 * t271) * r_i_i_C(1) + (t243 * t244 + t247 * t271) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:03
	% EndTime: 2019-10-10 11:07:03
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (378->75), mult. (797->120), div. (0->0), fcn. (746->10), ass. (0->56)
	t285 = sin(qJ(4));
	t326 = t285 * pkin(4);
	t288 = cos(qJ(4));
	t325 = t288 * pkin(4) + pkin(3) + pkin(8);
	t284 = cos(pkin(6));
	t289 = cos(qJ(2));
	t290 = cos(qJ(1));
	t317 = t290 * t289;
	t306 = t284 * t317;
	t286 = sin(qJ(2));
	t287 = sin(qJ(1));
	t320 = t287 * t286;
	t269 = -t306 + t320;
	t281 = qJD(4) + qJD(5);
	t324 = t269 * t281;
	t283 = sin(pkin(6));
	t323 = t283 * t287;
	t322 = t283 * t289;
	t321 = t283 * t290;
	t319 = t287 * t289;
	t318 = t290 * t286;
	t282 = qJ(4) + qJ(5);
	t279 = sin(t282);
	t280 = cos(t282);
	t313 = qJD(1) * t287;
	t305 = t283 * t313;
	t296 = -t305 - t324;
	t271 = t284 * t319 + t318;
	t297 = t284 * t318 + t319;
	t267 = t271 * qJD(1) + t297 * qJD(2);
	t298 = t281 * t321 + t267;
	t316 = (t296 * t279 + t298 * t280) * r_i_i_C(1) + (-t298 * t279 + t296 * t280) * r_i_i_C(2);
	t312 = qJD(1) * t290;
	t304 = t283 * t312;
	t295 = t271 * t281 + t304;
	t302 = qJD(2) * t284 + qJD(1);
	t265 = -qJD(1) * t306 - qJD(2) * t317 + t302 * t320;
	t299 = -t281 * t323 - t265;
	t261 = -t295 * t279 + t299 * t280;
	t262 = t299 * t279 + t295 * t280;
	t315 = t261 * r_i_i_C(1) - t262 * r_i_i_C(2);
	t311 = qJD(2) * t286;
	t303 = t283 * t311;
	t294 = -t281 * t284 + t303;
	t308 = t281 * t322;
	t314 = (t279 * t308 + t294 * t280) * r_i_i_C(1) + (-t294 * t279 + t280 * t308) * r_i_i_C(2);
	t310 = qJD(4) * t288;
	t309 = -r_i_i_C(3) - pkin(10) - pkin(9) - pkin(2);
	t307 = t284 * t320;
	t301 = r_i_i_C(1) * t280 - r_i_i_C(2) * t279;
	t300 = -t279 * r_i_i_C(1) - t280 * r_i_i_C(2);
	t293 = qJ(3) - t300 + t326;
	t292 = pkin(4) * t310 + t301 * t281 + qJD(3);
	t268 = -qJD(1) * t307 - t287 * t311 + t302 * t317;
	t266 = t297 * qJD(1) + t271 * qJD(2);
	t1 = [(-t267 * t279 - t280 * t324) * r_i_i_C(1) + (-t267 * t280 + t279 * t324) * r_i_i_C(2) - t267 * qJ(3) - t269 * qJD(3) - pkin(1) * t312 + (-t267 * t285 - t269 * t310) * pkin(4) + t309 * t268 + ((-qJD(4) * t326 + t300 * t281) * t290 + (-t301 - t325) * t313) * t283, -t309 * t265 + t292 * (-t307 + t317) - t293 * t266, -t265, (-t285 * t304 - t265 * t288 + (-t271 * t285 - t288 * t323) * qJD(4)) * pkin(4) + t315, t315, 0; t262 * r_i_i_C(1) + t261 * r_i_i_C(2) - t265 * qJ(3) + t271 * qJD(3) + t309 * t266 + (-pkin(1) * t287 + t325 * t321) * qJD(1) + (-t265 * t285 + (t271 * t288 - t285 * t323) * qJD(4)) * pkin(4), t309 * t267 + t293 * t268 + t292 * t297, t267, (-t285 * t305 + t267 * t288 + (-t269 * t285 + t288 * t321) * qJD(4)) * pkin(4) + t316, t316, 0; 0, (t292 * t286 + (t309 * t286 + t293 * t289) * qJD(2)) * t283, t303, (t288 * t303 + (-t284 * t288 + t285 * t322) * qJD(4)) * pkin(4) + t314, t314, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:05
	% EndTime: 2019-10-10 11:07:07
	% DurationCPUTime: 1.20s
	% Computational Cost: add. (999->116), mult. (1906->184), div. (0->0), fcn. (1888->12), ass. (0->75)
	t445 = sin(qJ(6));
	t449 = cos(qJ(6));
	t508 = t445 * r_i_i_C(1) + t449 * r_i_i_C(2);
	t506 = t508 * qJD(6);
	t501 = pkin(11) + r_i_i_C(3);
	t467 = t449 * r_i_i_C(1) - t445 * r_i_i_C(2);
	t464 = pkin(5) + t467;
	t444 = cos(pkin(6));
	t451 = cos(qJ(2));
	t452 = cos(qJ(1));
	t486 = t452 * t451;
	t475 = t444 * t486;
	t447 = sin(qJ(2));
	t448 = sin(qJ(1));
	t489 = t448 * t447;
	t426 = -t475 + t489;
	t442 = qJ(4) + qJ(5);
	t439 = sin(t442);
	t440 = cos(t442);
	t443 = sin(pkin(6));
	t490 = t443 * t452;
	t476 = t440 * t490;
	t422 = -t426 * t439 + t476;
	t487 = t452 * t447;
	t488 = t448 * t451;
	t427 = t444 * t487 + t488;
	t505 = -t422 * t445 - t427 * t449;
	t504 = t422 * t449 - t427 * t445;
	t462 = t467 * qJD(6);
	t428 = t444 * t488 + t487;
	t417 = qJD(1) * t428 + qJD(2) * t427;
	t441 = qJD(4) + qJD(5);
	t485 = qJD(1) * t443;
	t474 = t448 * t485;
	t461 = t426 * t441 + t474;
	t502 = -t417 * t440 + t439 * t461 - t441 * t476;
	t405 = t439 * (t441 * t490 + t417) + t440 * t461;
	t446 = sin(qJ(4));
	t458 = t446 * pkin(4) + t464 * t439 - t440 * t501 + qJ(3);
	t450 = cos(qJ(4));
	t498 = t450 * pkin(4);
	t497 = -pkin(2) - pkin(10) - pkin(9);
	t496 = pkin(8) + pkin(3) + t498;
	t493 = t443 * t447;
	t492 = t443 * t448;
	t491 = t443 * t451;
	t484 = qJD(2) * t447;
	t483 = qJD(2) * t451;
	t480 = t439 * t491;
	t479 = t439 * t492;
	t478 = t440 * t492;
	t477 = t444 * t489;
	t473 = t452 * t485;
	t472 = t443 * t483;
	t471 = t443 * t484;
	t468 = qJD(2) * t444 + qJD(1);
	t463 = t444 * t439 + t440 * t491;
	t460 = t428 * t441 + t473;
	t459 = t508 - t497;
	t415 = -qJD(1) * t475 - t452 * t483 + t468 * t489;
	t407 = t415 * t440 + t439 * t460 + t441 * t478;
	t408 = -t415 * t439 + t440 * t460 - t441 * t479;
	t457 = -t506 * (t428 * t440 - t479) + t501 * t408 - t464 * t407;
	t456 = -t506 * (t426 * t440 + t439 * t490) + t501 * t405 - t464 * t502;
	t413 = -t439 * t471 + t441 * t463;
	t455 = t506 * t463 + t464 * (t441 * t480 + (-t441 * t444 + t471) * t440) - t501 * t413;
	t454 = qJD(4) * t498 + qJD(3) - t439 * t506 + (t439 * t501 + t464 * t440) * t441;
	t429 = -t477 + t486;
	t425 = t444 * t440 - t480;
	t420 = t428 * t439 + t478;
	t418 = -qJD(1) * t477 - t448 * t484 + t468 * t486;
	t416 = qJD(1) * t427 + qJD(2) * t428;
	t396 = t408 * t449 - t416 * t445 + (-t420 * t445 + t429 * t449) * qJD(6);
	t395 = -t408 * t445 - t416 * t449 + (-t420 * t449 - t429 * t445) * qJD(6);
	t1 = [-t417 * qJ(3) - t426 * qJD(3) - t464 * t405 - t459 * t418 - t501 * t502 + (t505 * r_i_i_C(1) - t504 * r_i_i_C(2)) * qJD(6) + (-t452 * pkin(1) - t492 * t496) * qJD(1) + (-t417 * t446 + (-t426 * t450 - t446 * t490) * qJD(4)) * pkin(4), t459 * t415 - t458 * t416 - t428 * t462 + t454 * t429, -t415, (-t446 * t473 - t415 * t450 + (-t428 * t446 - t450 * t492) * qJD(4)) * pkin(4) + t457, t457, t395 * r_i_i_C(1) - t396 * r_i_i_C(2); t408 * pkin(5) + t396 * r_i_i_C(1) + t395 * r_i_i_C(2) - t415 * qJ(3) + t428 * qJD(3) + t497 * t416 + t501 * t407 + (-pkin(1) * t448 + t490 * t496) * qJD(1) + (-t415 * t446 + (t428 * t450 - t446 * t492) * qJD(4)) * pkin(4), -t417 * t459 + t418 * t458 - t426 * t462 + t427 * t454, t417, (-t446 * t474 + t417 * t450 + (-t426 * t446 + t450 * t490) * qJD(4)) * pkin(4) + t456, t456, (-t405 * t445 + t418 * t449) * r_i_i_C(1) + (-t405 * t449 - t418 * t445) * r_i_i_C(2) + (t504 * r_i_i_C(1) + t505 * r_i_i_C(2)) * qJD(6); 0, ((qJD(2) * t458 + t462) * t451 + (-qJD(2) * t459 + t454) * t447) * t443, t471, (t450 * t471 + (-t444 * t450 + t446 * t491) * qJD(4)) * pkin(4) + t455, t455, (t413 * t445 + t449 * t472) * r_i_i_C(1) + (t413 * t449 - t445 * t472) * r_i_i_C(2) + ((-t425 * t449 - t445 * t493) * r_i_i_C(1) + (t425 * t445 - t449 * t493) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end