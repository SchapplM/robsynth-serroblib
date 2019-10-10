% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6PPRRRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:22
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRRR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobiaD_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:21
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (151->10), mult. (481->28), div. (18->4), fcn. (595->10), ass. (0->20)
	t84 = sin(pkin(14));
	t88 = cos(pkin(14));
	t89 = cos(pkin(13));
	t85 = sin(pkin(13));
	t98 = t85 * cos(pkin(6));
	t82 = -t84 * t98 + t89 * t88;
	t92 = sin(qJ(3));
	t93 = cos(qJ(3));
	t96 = t85 * sin(pkin(7)) * sin(pkin(6)) + (-t89 * t84 - t88 * t98) * cos(pkin(7));
	t78 = t82 * t92 - t96 * t93;
	t79 = t82 * t93 + t96 * t92;
	t76 = 0.1e1 / t79 ^ 2;
	t104 = qJD(3) * t76;
	t101 = t79 * t104;
	t102 = t78 / t79 * t104;
	t75 = t78 ^ 2;
	t72 = t75 * t76 + 0.1e1;
	t103 = (t78 * t101 + t75 * t102) / t72 ^ 2;
	t70 = 0.1e1 / t72;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, -0.2e1 * t103 + 0.2e1 * (t70 * t101 + (t70 * t102 - t76 * t103) * t78) * t78, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:21
	% EndTime: 2019-10-09 21:22:22
	% DurationCPUTime: 1.33s
	% Computational Cost: add. (4180->76), mult. (12611->185), div. (275->12), fcn. (16541->17), ass. (0->91)
	t235 = sin(pkin(14));
	t240 = cos(pkin(14));
	t241 = cos(pkin(13));
	t236 = sin(pkin(13));
	t244 = cos(pkin(6));
	t272 = t236 * t244;
	t232 = -t235 * t272 + t241 * t240;
	t246 = sin(qJ(3));
	t248 = cos(qJ(3));
	t231 = -t241 * t235 - t240 * t272;
	t243 = cos(pkin(7));
	t238 = sin(pkin(7));
	t239 = sin(pkin(6));
	t271 = t238 * t239;
	t258 = t231 * t243 + t236 * t271;
	t292 = t232 * t246 - t258 * t248;
	t268 = t240 * t243;
	t270 = t238 * t244;
	t290 = (-t235 * t246 + t248 * t268) * t239 + t248 * t270;
	t221 = t232 * t248 + t258 * t246;
	t267 = t241 * t244;
	t230 = t235 * t267 + t236 * t240;
	t257 = -t236 * t235 + t240 * t267;
	t255 = -t241 * t271 + t257 * t243;
	t289 = -t230 * t246 + t255 * t248;
	t237 = sin(pkin(8));
	t242 = cos(pkin(8));
	t269 = t239 * t243;
	t209 = t289 * t237 - (-t257 * t238 - t241 * t269) * t242;
	t222 = -t290 * t237 + (-t240 * t271 + t244 * t243) * t242;
	t204 = atan2(t209, t222);
	t199 = sin(t204);
	t200 = cos(t204);
	t186 = t199 * t209 + t200 * t222;
	t183 = 0.1e1 / t186;
	t245 = sin(qJ(4));
	t247 = cos(qJ(4));
	t228 = -t231 * t238 + t236 * t269;
	t277 = t228 * t237;
	t260 = -t242 * t292 + t277;
	t198 = t221 * t247 + t260 * t245;
	t194 = 0.1e1 / t198;
	t216 = 0.1e1 / t222;
	t234 = t237 ^ 2;
	t184 = 0.1e1 / t186 ^ 2;
	t195 = 0.1e1 / t198 ^ 2;
	t217 = 0.1e1 / t222 ^ 2;
	t210 = t228 * t242 + t237 * t292;
	t208 = t210 ^ 2;
	t182 = t208 * t184 + 0.1e1;
	t215 = t221 * qJD(3);
	t219 = -t230 * t248 - t255 * t246;
	t213 = t219 * qJD(3);
	t227 = -t246 * t270 + (-t235 * t248 - t246 * t268) * t239;
	t225 = t227 * qJD(3);
	t281 = t209 * t217;
	t207 = t209 ^ 2;
	t203 = t207 * t217 + 0.1e1;
	t201 = 0.1e1 / t203;
	t282 = t201 * t237;
	t178 = (t213 * t216 + t225 * t281) * t282;
	t261 = -t199 * t222 + t200 * t209;
	t175 = (t199 * t213 - t200 * t225) * t237 + t261 * t178;
	t287 = t175 * t183 * t184;
	t288 = (t210 * t184 * t215 * t237 - t208 * t287) / t182 ^ 2;
	t265 = t242 * t247;
	t278 = t221 * t245;
	t197 = -t247 * t277 + t265 * t292 + t278;
	t193 = t197 ^ 2;
	t190 = t193 * t195 + 0.1e1;
	t214 = t292 * qJD(3);
	t266 = t242 * t245;
	t192 = -t215 * t266 - t214 * t247 + (t260 * t247 - t278) * qJD(4);
	t284 = t192 * t194 * t195;
	t191 = t198 * qJD(4) - t214 * t245 + t215 * t265;
	t285 = t191 * t195;
	t286 = (-t193 * t284 + t197 * t285) / t190 ^ 2;
	t206 = -t221 * t266 - t247 * t292;
	t283 = t197 * t206;
	t280 = t209 * t227;
	t279 = t216 * t217 * t225;
	t259 = t216 * t219 + t217 * t280;
	t205 = t221 * t265 - t245 * t292;
	t224 = t290 * qJD(3);
	t212 = t289 * qJD(3);
	t188 = 0.1e1 / t190;
	t180 = 0.1e1 / t182;
	t179 = t259 * t282;
	t176 = (t199 * t219 - t200 * t227) * t237 + t261 * t179;
	t174 = -0.2e1 * t259 * t234 / t203 ^ 2 * (t207 * t279 + t213 * t281) + (0.2e1 * t234 * t279 * t280 - t212 * t216 * t237 + (-t209 * t224 * t237 + (t213 * t227 + t219 * t225) * t234) * t217) * t201;
	t1 = [0, 0, t174, 0, 0, 0; 0, 0, (-(-t186 * t179 * t178 + t261 * t174) * t184 * t180 + 0.2e1 * (t180 * t287 + t184 * t288) * t176) * t210 + (-0.2e1 * t221 * t183 * t288 + (-t214 * t183 + (-t221 * t175 - t176 * t215 + (-(t178 * t219 + t179 * t213 + t224) * t200 - (t178 * t227 + t179 * t225 - t212) * t199) * t210) * t184) * t180) * t237, 0, 0, 0; 0, 0, 0.2e1 * (-t194 * t205 + t195 * t283) * t286 + ((qJD(4) * t206 - t214 * t265 - t215 * t245) * t194 + 0.2e1 * t283 * t284 + (-t205 * t192 - (-qJD(4) * t205 + t214 * t266 - t215 * t247) * t197 - t206 * t191) * t195) * t188, -0.2e1 * t286 + 0.2e1 * (t188 * t285 + (-t188 * t284 - t195 * t286) * t197) * t197, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:22
	% EndTime: 2019-10-09 21:22:25
	% DurationCPUTime: 3.67s
	% Computational Cost: add. (15480->143), mult. (47073->274), div. (538->12), fcn. (61235->19), ass. (0->130)
	t446 = sin(pkin(14));
	t447 = sin(pkin(13));
	t411 = t447 * t446;
	t450 = cos(pkin(14));
	t451 = cos(pkin(13));
	t417 = t451 * t450;
	t453 = cos(pkin(6));
	t391 = -t453 * t417 + t411;
	t448 = sin(pkin(7));
	t449 = sin(pkin(6));
	t415 = t449 * t448;
	t452 = cos(pkin(7));
	t460 = -t391 * t452 - t451 * t415;
	t413 = t447 * t450;
	t416 = t451 * t446;
	t360 = -t453 * t413 - t416;
	t412 = t447 * t449;
	t459 = t452 * t360 + t448 * t412;
	t418 = t452 * t449;
	t458 = t450 * t418 + t453 * t448;
	t368 = sin(qJ(4));
	t369 = sin(qJ(3));
	t392 = t453 * t416 + t413;
	t455 = cos(qJ(3));
	t381 = -t392 * t369 + t460 * t455;
	t454 = cos(qJ(4));
	t378 = t381 * t454;
	t348 = t460 * t369 + t392 * t455;
	t380 = qJD(3) * t348;
	t457 = qJD(4) * t378 - t368 * t380;
	t366 = cos(pkin(8));
	t365 = sin(pkin(8));
	t383 = (t391 * t448 - t451 * t418) * t365;
	t325 = t348 * t454 + (t381 * t366 + t383) * t368;
	t345 = t381 * qJD(3);
	t425 = t366 * t454;
	t299 = t325 * qJD(4) + t345 * t368 + t380 * t425;
	t382 = t454 * t383;
	t436 = t348 * t368;
	t323 = -t366 * t378 - t382 + t436;
	t321 = t323 ^ 2;
	t414 = t449 * t446;
	t356 = t458 * t369 + t455 * t414;
	t355 = -t369 * t414 + t458 * t455;
	t424 = t454 * t355;
	t435 = (-t450 * t415 + t453 * t452) * t365;
	t393 = -t356 * t368 + t366 * t424 + t454 * t435;
	t334 = 0.1e1 / t393 ^ 2;
	t313 = t321 * t334 + 0.1e1;
	t437 = t323 * t334;
	t337 = t356 * t454 + (t355 * t366 + t435) * t368;
	t353 = t355 * qJD(3);
	t354 = t356 * qJD(3);
	t319 = t337 * qJD(4) + t353 * t368 + t354 * t425;
	t333 = 0.1e1 / t393;
	t438 = t319 * t333 * t334;
	t456 = -0.2e1 * (t299 * t437 + t321 * t438) / t313 ^ 2;
	t361 = -t453 * t411 + t417;
	t349 = -t361 * t369 + t459 * t455;
	t350 = t361 * t455 + t459 * t369;
	t314 = atan2(-t323, -t393);
	t309 = sin(t314);
	t310 = cos(t314);
	t293 = -t309 * t323 - t310 * t393;
	t290 = 0.1e1 / t293;
	t395 = -t360 * t448 + t452 * t412;
	t389 = t395 * t365;
	t327 = t350 * t454 + (t349 * t366 + t389) * t368;
	t338 = -t349 * t365 + t395 * t366;
	t367 = sin(qJ(5));
	t370 = cos(qJ(5));
	t308 = t327 * t370 + t338 * t367;
	t304 = 0.1e1 / t308;
	t291 = 0.1e1 / t293 ^ 2;
	t305 = 0.1e1 / t308 ^ 2;
	t326 = -t349 * t425 + t350 * t368 - t454 * t389;
	t322 = t326 ^ 2;
	t289 = t291 * t322 + 0.1e1;
	t346 = t349 * qJD(3);
	t347 = t350 * qJD(3);
	t301 = t327 * qJD(4) + t346 * t368 + t347 * t425;
	t441 = t291 * t326;
	t311 = 0.1e1 / t313;
	t283 = (t299 * t333 + t319 * t437) * t311;
	t410 = t309 * t393 - t310 * t323;
	t279 = t410 * t283 - t309 * t299 + t310 * t319;
	t444 = t279 * t290 * t291;
	t445 = (t301 * t441 - t322 * t444) / t289 ^ 2;
	t431 = t366 * t368;
	t302 = -t326 * qJD(4) + t346 * t454 - t347 * t431;
	t432 = t365 * t370;
	t294 = t308 * qJD(5) + t302 * t367 - t347 * t432;
	t307 = t327 * t367 - t338 * t370;
	t303 = t307 ^ 2;
	t298 = t303 * t305 + 0.1e1;
	t439 = t305 * t307;
	t430 = qJD(5) * t307;
	t433 = t365 * t367;
	t295 = t302 * t370 + t347 * t433 - t430;
	t440 = t295 * t304 * t305;
	t443 = (t294 * t439 - t303 * t440) / t298 ^ 2;
	t287 = 0.1e1 / t289;
	t442 = t287 * t291;
	t429 = 0.2e1 * t445;
	t428 = -0.2e1 * t443;
	t427 = t307 * t440;
	t423 = qJD(4) * t436;
	t421 = 0.2e1 * t326 * t444;
	t420 = 0.2e1 * t323 * t438;
	t407 = -t304 * t367 + t370 * t439;
	t406 = t325 * t333 + t337 * t437;
	t329 = t348 * t425 + t381 * t368;
	t339 = t355 * t368 + t356 * t425;
	t405 = t329 * t333 + t339 * t437;
	t331 = t349 * t454 - t350 * t431;
	t404 = -t331 * t367 + t350 * t432;
	t318 = t331 * t370 + t350 * t433;
	t398 = -t349 * t368 - t350 * t425;
	t328 = t353 * t425 - t354 * t368 + (-t356 * t431 + t424) * qJD(4);
	t320 = t393 * qJD(4) + t353 * t454 - t354 * t431;
	t316 = t398 * qJD(4) - t346 * t431 - t347 * t454;
	t315 = t345 * t425 - t366 * t423 + t457;
	t300 = qJD(4) * t382 + t345 * t454 + t457 * t366 - t423;
	t296 = 0.1e1 / t298;
	t285 = t405 * t311;
	t284 = t406 * t311;
	t280 = t410 * t284 - t309 * t325 + t310 * t337;
	t278 = t405 * t456 + (t339 * t420 + t315 * t333 + (t299 * t339 + t319 * t329 + t323 * t328) * t334) * t311;
	t276 = t406 * t456 + (t337 * t420 + t300 * t333 + (t299 * t337 + t319 * t325 + t320 * t323) * t334) * t311;
	t1 = [0, 0, t278, t276, 0, 0; 0, 0, -(-t279 * t442 - 0.2e1 * t290 * t445) * t398 + ((t331 * qJD(4) + t346 * t425 - t347 * t368) * t290 - ((-t278 * t323 - t285 * t299 + t328 + (t285 * t393 - t329) * t283) * t310 + (t278 * t393 - t285 * t319 - t315 + (t285 * t323 - t339) * t283) * t309) * t441) * t287 + (t287 * t421 - t301 * t442 + t441 * t429) * (t410 * t285 - t309 * t329 + t310 * t339), (t280 * t441 - t290 * t327) * t429 + (t280 * t421 + t302 * t290 + (-t327 * t279 - t280 * t301 + (-(-t276 * t323 - t284 * t299 + t320 + (t284 * t393 - t325) * t283) * t310 - (t276 * t393 - t284 * t319 - t300 + (t284 * t323 - t337) * t283) * t309) * t326) * t291) * t287, 0, 0; 0, 0, 0.2e1 * (t304 * t404 + t318 * t439) * t443 + ((t318 * qJD(5) + t316 * t367 - t346 * t432) * t304 + 0.2e1 * t318 * t427 + (t404 * t295 - (t404 * qJD(5) + t316 * t370 + t346 * t433) * t307 - t318 * t294) * t305) * t296, t407 * t326 * t428 + (t407 * t301 + ((-qJD(5) * t304 - 0.2e1 * t427) * t370 + (t294 * t370 + (t295 - t430) * t367) * t305) * t326) * t296, t428 + 0.2e1 * (t294 * t296 * t305 + (-t296 * t440 - t305 * t443) * t307) * t307, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:23
	% EndTime: 2019-10-09 21:22:31
	% DurationCPUTime: 8.36s
	% Computational Cost: add. (40633->214), mult. (120610->400), div. (816->12), fcn. (158008->21), ass. (0->170)
	t574 = sin(pkin(14));
	t575 = sin(pkin(13));
	t533 = t575 * t574;
	t578 = cos(pkin(14));
	t579 = cos(pkin(13));
	t539 = t579 * t578;
	t581 = cos(pkin(6));
	t513 = -t581 * t539 + t533;
	t576 = sin(pkin(7));
	t577 = sin(pkin(6));
	t537 = t577 * t576;
	t580 = cos(pkin(7));
	t590 = t513 * t580 + t579 * t537;
	t535 = t575 * t578;
	t538 = t579 * t574;
	t515 = t581 * t535 + t538;
	t534 = t575 * t577;
	t589 = t515 * t580 - t576 * t534;
	t540 = t580 * t577;
	t588 = t578 * t540 + t581 * t576;
	t482 = sin(qJ(4));
	t486 = cos(qJ(4));
	t478 = sin(pkin(8));
	t479 = cos(pkin(8));
	t483 = sin(qJ(3));
	t514 = t581 * t538 + t535;
	t582 = cos(qJ(3));
	t498 = t514 * t483 + t590 * t582;
	t502 = t513 * t576 - t579 * t540;
	t493 = t502 * t478 - t498 * t479;
	t586 = t590 * t483 - t514 * t582;
	t442 = t482 * t586 + t493 * t486;
	t459 = t498 * qJD(3);
	t460 = t586 * qJD(3);
	t556 = t479 * t482;
	t421 = t442 * qJD(4) - t459 * t486 + t460 * t556;
	t443 = t493 * t482 - t486 * t586;
	t481 = sin(qJ(5));
	t485 = cos(qJ(5));
	t494 = t498 * t478 + t502 * t479;
	t428 = t443 * t485 + t494 * t481;
	t557 = t478 * t485;
	t395 = t428 * qJD(5) + t421 * t481 + t460 * t557;
	t426 = t443 * t481 - t494 * t485;
	t424 = t426 ^ 2;
	t536 = t577 * t574;
	t470 = t588 * t483 + t582 * t536;
	t469 = -t483 * t536 + t588 * t582;
	t473 = -t578 * t537 + t581 * t580;
	t531 = t469 * t479 + t473 * t478;
	t455 = t470 * t486 + t531 * t482;
	t465 = -t469 * t478 + t473 * t479;
	t446 = t455 * t481 - t465 * t485;
	t440 = 0.1e1 / t446 ^ 2;
	t411 = t424 * t440 + 0.1e1;
	t409 = 0.1e1 / t411;
	t454 = -t470 * t482 + t531 * t486;
	t467 = t469 * qJD(3);
	t468 = t470 * qJD(3);
	t437 = t454 * qJD(4) + t467 * t486 - t468 * t556;
	t447 = t455 * t485 + t465 * t481;
	t413 = t447 * qJD(5) + t437 * t481 - t468 * t557;
	t439 = 0.1e1 / t446;
	t564 = t426 * t440;
	t378 = (-t395 * t439 + t413 * t564) * t409;
	t412 = atan2(-t426, t446);
	t407 = sin(t412);
	t408 = cos(t412);
	t532 = -t407 * t446 - t408 * t426;
	t373 = t532 * t378 - t407 * t395 + t408 * t413;
	t391 = -t407 * t426 + t408 * t446;
	t388 = 0.1e1 / t391;
	t389 = 0.1e1 / t391 ^ 2;
	t587 = t373 * t388 * t389;
	t516 = -t581 * t533 + t539;
	t464 = -t589 * t483 + t516 * t582;
	t499 = t516 * t483 + t589 * t582;
	t503 = t515 * t576 + t580 * t534;
	t500 = t503 * t478;
	t495 = -t499 * t479 + t500;
	t445 = t464 * t486 + t495 * t482;
	t496 = t499 * t478 + t503 * t479;
	t429 = t445 * t481 - t496 * t485;
	t545 = 0.2e1 * t429 * t587;
	t526 = -t439 * t442 + t454 * t564;
	t585 = t481 * t526;
	t565 = t413 * t439 * t440;
	t583 = -0.2e1 * (t395 * t564 - t424 * t565) / t411 ^ 2;
	t430 = t445 * t485 + t496 * t481;
	t497 = t499 * t486;
	t560 = t464 * t482;
	t444 = t479 * t497 - t486 * t500 + t560;
	t480 = sin(qJ(6));
	t484 = cos(qJ(6));
	t406 = t430 * t484 + t444 * t480;
	t402 = 0.1e1 / t406;
	t403 = 0.1e1 / t406 ^ 2;
	t461 = t499 * qJD(3);
	t462 = t464 * qJD(3);
	t423 = -t462 * t556 - t461 * t486 + (t495 * t486 - t560) * qJD(4);
	t558 = t478 * t481;
	t398 = -t429 * qJD(5) + t423 * t485 + t462 * t558;
	t555 = t479 * t486;
	t422 = t445 * qJD(4) - t461 * t482 + t462 * t555;
	t386 = t406 * qJD(6) + t398 * t480 - t422 * t484;
	t405 = t430 * t480 - t444 * t484;
	t401 = t405 ^ 2;
	t394 = t401 * t403 + 0.1e1;
	t569 = t403 * t405;
	t551 = qJD(6) * t405;
	t387 = t398 * t484 + t422 * t480 - t551;
	t572 = t387 * t402 * t403;
	t573 = (t386 * t569 - t401 * t572) / t394 ^ 2;
	t571 = t389 * t429;
	t397 = t430 * qJD(5) + t423 * t481 - t462 * t557;
	t570 = t397 * t389;
	t568 = t405 * t484;
	t567 = t407 * t429;
	t566 = t408 * t429;
	t563 = t444 * t481;
	t562 = t444 * t485;
	t554 = t480 * t402;
	t553 = qJD(5) * t481;
	t552 = qJD(5) * t485;
	t425 = t429 ^ 2;
	t385 = t425 * t389 + 0.1e1;
	t550 = 0.2e1 * (-t425 * t587 + t429 * t570) / t385 ^ 2;
	t549 = -0.2e1 * t573;
	t548 = 0.2e1 * t573;
	t546 = t405 * t572;
	t544 = 0.2e1 * t546;
	t543 = -0.2e1 * t426 * t565;
	t542 = qJD(6) * t562 + t423;
	t450 = -t464 * t556 - t497;
	t435 = t450 * t485 + t464 * t558;
	t449 = t464 * t555 - t499 * t482;
	t418 = t435 * t484 + t449 * t480;
	t417 = t435 * t480 - t449 * t484;
	t529 = t403 * t568 - t554;
	t528 = -t428 * t439 + t447 * t564;
	t448 = -t498 * t486 + t556 * t586;
	t433 = t448 * t481 + t557 * t586;
	t456 = t469 * t486 - t470 * t556;
	t451 = t456 * t481 - t470 * t557;
	t527 = -t433 * t439 + t451 * t564;
	t525 = -t450 * t481 + t464 * t557;
	t519 = qJD(6) * t445 - t422 * t485 + t444 * t553;
	t436 = -t455 * qJD(4) - t467 * t482 - t468 * t555;
	t432 = -t449 * qJD(4) + t461 * t556 - t462 * t486;
	t431 = t450 * qJD(4) - t461 * t555 - t462 * t482;
	t420 = -t443 * qJD(4) + t459 * t482 + t460 * t555;
	t419 = (-t467 * t556 - t468 * t486 + (-t469 * t482 - t470 * t555) * qJD(4)) * t481 - t467 * t557 + (t456 * t485 + t470 * t558) * qJD(5);
	t416 = t445 * t480 - t484 * t562;
	t415 = -t445 * t484 - t480 * t562;
	t414 = -t446 * qJD(5) + t437 * t485 + t468 * t558;
	t400 = t525 * qJD(5) + t432 * t485 - t461 * t558;
	t399 = (t459 * t556 + t460 * t486 + (t498 * t482 + t555 * t586) * qJD(4)) * t481 + t448 * t552 + t459 * t557 - t586 * t478 * t553;
	t396 = -t426 * qJD(5) + t421 * t485 - t460 * t558;
	t392 = 0.1e1 / t394;
	t383 = 0.1e1 / t385;
	t382 = t409 * t585;
	t381 = t527 * t409;
	t380 = t528 * t409;
	t376 = (-t407 * t442 + t408 * t454) * t481 + t532 * t382;
	t375 = t532 * t381 - t407 * t433 + t408 * t451;
	t374 = t532 * t380 - t407 * t428 + t408 * t447;
	t372 = t527 * t583 + (t451 * t543 - t399 * t439 + (t395 * t451 + t413 * t433 + t419 * t426) * t440) * t409;
	t370 = t528 * t583 + (t447 * t543 - t396 * t439 + (t395 * t447 + t413 * t428 + t414 * t426) * t440) * t409;
	t369 = t583 * t585 + (t526 * t552 + (t454 * t543 - t420 * t439 + (t395 * t454 + t413 * t442 + t426 * t436) * t440) * t481) * t409;
	t1 = [0, 0, t372, t369, t370, 0; 0, 0, (t375 * t571 + t388 * t525) * t550 + ((t435 * qJD(5) + t432 * t481 + t461 * t557) * t388 + t375 * t545 + (t525 * t373 - t375 * t397 - (-t372 * t426 - t381 * t395 + t419 + (-t381 * t446 - t433) * t378) * t566 - (-t372 * t446 - t381 * t413 - t399 + (t381 * t426 - t451) * t378) * t567) * t389) * t383, (t376 * t571 + t388 * t563) * t550 + ((-t422 * t481 - t444 * t552) * t388 + (-t570 + t545) * t376 + (t563 * t373 - (t454 * t552 - t369 * t426 - t382 * t395 + t436 * t481 + (-t382 * t446 - t442 * t481) * t378) * t566 - (-t442 * t552 - t369 * t446 - t382 * t413 - t420 * t481 + (t382 * t426 - t454 * t481) * t378) * t567) * t389) * t383, (t374 * t571 - t388 * t430) * t550 + (t374 * t545 + t398 * t388 + (-t430 * t373 - t374 * t397 - (-t370 * t426 - t380 * t395 + t414 + (-t380 * t446 - t428) * t378) * t566 - (-t370 * t446 - t380 * t413 - t396 + (t380 * t426 - t447) * t378) * t567) * t389) * t383, 0; 0, 0, (-t402 * t417 + t418 * t569) * t548 + ((t418 * qJD(6) + t400 * t480 - t431 * t484) * t402 + t418 * t544 + (-t417 * t387 - (-t417 * qJD(6) + t400 * t484 + t431 * t480) * t405 - t418 * t386) * t403) * t392, (-t402 * t415 + t416 * t569) * t548 + (t416 * t544 - t542 * t402 * t484 + t519 * t554 + (-t542 * t405 * t480 - t416 * t386 - t415 * t387 - t519 * t568) * t403) * t392, t529 * t429 * t549 + (t529 * t397 + ((-qJD(6) * t402 - 0.2e1 * t546) * t484 + (t386 * t484 + (t387 - t551) * t480) * t403) * t429) * t392, t549 + 0.2e1 * (t386 * t403 * t392 + (-t392 * t572 - t403 * t573) * t405) * t405;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end