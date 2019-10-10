% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRRP2
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
%   Wie in S6PPRRRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:16
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRRP2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRP2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:16:21
	% EndTime: 2019-10-09 21:16:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:16:21
	% EndTime: 2019-10-09 21:16:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:16:21
	% EndTime: 2019-10-09 21:16:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:16:21
	% EndTime: 2019-10-09 21:16:21
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (151->10), mult. (481->28), div. (18->4), fcn. (595->10), ass. (0->20)
	t84 = sin(pkin(12));
	t88 = cos(pkin(12));
	t89 = cos(pkin(11));
	t85 = sin(pkin(11));
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
	% StartTime: 2019-10-09 21:16:21
	% EndTime: 2019-10-09 21:16:22
	% DurationCPUTime: 0.97s
	% Computational Cost: add. (2591->65), mult. (8150->143), div. (281->12), fcn. (10563->15), ass. (0->78)
	t208 = sin(pkin(7));
	t209 = sin(pkin(6));
	t211 = cos(pkin(6));
	t207 = sin(pkin(12));
	t250 = sin(pkin(11));
	t234 = t250 * t207;
	t210 = cos(pkin(12));
	t251 = cos(pkin(11));
	t235 = t251 * t210;
	t252 = cos(pkin(7));
	t256 = (t211 * t235 - t234) * t252 - t208 * t209 * t251;
	t233 = t250 * t210;
	t236 = t251 * t207;
	t224 = t211 * t233 + t236;
	t238 = t209 * t250;
	t255 = -t208 * t238 + t224 * t252;
	t202 = t211 * t236 + t233;
	t213 = sin(qJ(3));
	t253 = cos(qJ(3));
	t188 = t202 * t253 + t213 * t256;
	t237 = t210 * t252;
	t241 = t208 * t211;
	t254 = (-t207 * t213 + t253 * t237) * t209 + t253 * t241;
	t186 = t202 * t213 - t253 * t256;
	t179 = atan2(-t186, -t254);
	t174 = sin(t179);
	t175 = cos(t179);
	t162 = -t174 * t186 - t175 * t254;
	t159 = 0.1e1 / t162;
	t203 = -t211 * t234 + t235;
	t190 = t203 * t253 - t255 * t213;
	t199 = t224 * t208 + t252 * t238;
	t212 = sin(qJ(4));
	t214 = cos(qJ(4));
	t173 = t190 * t214 + t199 * t212;
	t169 = 0.1e1 / t173;
	t194 = 0.1e1 / t254;
	t160 = 0.1e1 / t162 ^ 2;
	t170 = 0.1e1 / t173 ^ 2;
	t195 = 0.1e1 / t254 ^ 2;
	t184 = t186 ^ 2;
	t178 = t184 * t195 + 0.1e1;
	t176 = 0.1e1 / t178;
	t181 = t188 * qJD(3);
	t198 = t213 * t241 + (t253 * t207 + t213 * t237) * t209;
	t192 = t198 * qJD(3);
	t244 = t186 * t195;
	t153 = (t181 * t194 + t192 * t244) * t176;
	t229 = t174 * t254 - t175 * t186;
	t150 = t229 * t153 - t174 * t181 + t175 * t192;
	t249 = t150 * t159 * t160;
	t172 = t190 * t212 - t199 * t214;
	t168 = t172 ^ 2;
	t165 = t168 * t170 + 0.1e1;
	t189 = t203 * t213 + t255 * t253;
	t182 = t189 * qJD(3);
	t166 = t173 * qJD(4) - t182 * t212;
	t245 = t170 * t172;
	t240 = qJD(4) * t172;
	t167 = -t182 * t214 - t240;
	t246 = t167 * t169 * t170;
	t248 = (t166 * t245 - t168 * t246) / t165 ^ 2;
	t247 = t160 * t189;
	t243 = t186 * t198;
	t242 = t192 * t194 * t195;
	t239 = -0.2e1 * t248;
	t227 = -t169 * t212 + t214 * t245;
	t226 = t188 * t194 + t195 * t243;
	t191 = t254 * qJD(3);
	t185 = t189 ^ 2;
	t183 = t190 * qJD(3);
	t180 = t186 * qJD(3);
	t163 = 0.1e1 / t165;
	t157 = t185 * t160 + 0.1e1;
	t154 = t226 * t176;
	t151 = t229 * t154 - t174 * t188 + t175 * t198;
	t149 = -0.2e1 * t226 / t178 ^ 2 * (t181 * t244 + t184 * t242) + (0.2e1 * t242 * t243 - t180 * t194 + (t181 * t198 + t186 * t191 + t188 * t192) * t195) * t176;
	t1 = [0, 0, t149, 0, 0, 0; 0, 0, 0.2e1 * (t151 * t247 - t159 * t190) / t157 ^ 2 * (t183 * t247 - t185 * t249) + (-t182 * t159 + (-t190 * t150 - t151 * t183) * t160 + (0.2e1 * t151 * t249 + (-(-t149 * t186 - t154 * t181 + t191 + (t154 * t254 - t188) * t153) * t175 - (t149 * t254 - t154 * t192 + t180 + (t154 * t186 - t198) * t153) * t174) * t160) * t189) / t157, 0, 0, 0; 0, 0, t227 * t189 * t239 + (t227 * t183 + ((-qJD(4) * t169 - 0.2e1 * t172 * t246) * t214 + (t166 * t214 + (t167 - t240) * t212) * t170) * t189) * t163, t239 + 0.2e1 * (t163 * t166 * t170 + (-t163 * t246 - t170 * t248) * t172) * t172, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:16:22
	% EndTime: 2019-10-09 21:16:24
	% DurationCPUTime: 2.50s
	% Computational Cost: add. (9293->119), mult. (27507->244), div. (559->12), fcn. (36133->17), ass. (0->119)
	t368 = sin(pkin(12));
	t369 = sin(pkin(11));
	t337 = t369 * t368;
	t372 = cos(pkin(12));
	t373 = cos(pkin(11));
	t343 = t373 * t372;
	t375 = cos(pkin(6));
	t308 = -t375 * t337 + t343;
	t312 = sin(qJ(3));
	t315 = cos(qJ(3));
	t339 = t369 * t372;
	t342 = t373 * t368;
	t329 = t375 * t339 + t342;
	t371 = sin(pkin(6));
	t338 = t369 * t371;
	t370 = sin(pkin(7));
	t374 = cos(pkin(7));
	t377 = t329 * t374 - t370 * t338;
	t292 = t308 * t312 + t377 * t315;
	t307 = t375 * t342 + t339;
	t328 = -t375 * t343 + t337;
	t341 = t371 * t370;
	t323 = -t328 * t374 - t373 * t341;
	t291 = t307 * t315 + t323 * t312;
	t311 = sin(qJ(4));
	t314 = cos(qJ(4));
	t344 = t374 * t371;
	t322 = t328 * t370 - t373 * t344;
	t282 = t291 * t314 + t322 * t311;
	t290 = -t307 * t312 + t323 * t315;
	t285 = t290 * qJD(3);
	t258 = t282 * qJD(4) + t285 * t311;
	t280 = t291 * t311 - t322 * t314;
	t278 = t280 ^ 2;
	t327 = t372 * t344 + t375 * t370;
	t340 = t371 * t368;
	t304 = t327 * t312 + t315 * t340;
	t306 = -t372 * t341 + t375 * t374;
	t297 = t304 * t311 - t306 * t314;
	t295 = 0.1e1 / t297 ^ 2;
	t272 = t278 * t295 + 0.1e1;
	t270 = 0.1e1 / t272;
	t298 = t304 * t314 + t306 * t311;
	t303 = -t312 * t340 + t327 * t315;
	t299 = t303 * qJD(3);
	t276 = t298 * qJD(4) + t299 * t311;
	t294 = 0.1e1 / t297;
	t357 = t280 * t295;
	t242 = (-t258 * t294 + t276 * t357) * t270;
	t273 = atan2(-t280, t297);
	t268 = sin(t273);
	t269 = cos(t273);
	t336 = -t268 * t297 - t269 * t280;
	t238 = t336 * t242 - t268 * t258 + t269 * t276;
	t252 = -t268 * t280 + t269 * t297;
	t249 = 0.1e1 / t252;
	t250 = 0.1e1 / t252 ^ 2;
	t380 = t238 * t249 * t250;
	t293 = t308 * t315 - t312 * t377;
	t324 = t329 * t370 + t374 * t338;
	t283 = t293 * t311 - t324 * t314;
	t379 = 0.2e1 * t283 * t380;
	t332 = -t290 * t294 + t303 * t357;
	t378 = t311 * t332;
	t358 = t276 * t294 * t295;
	t376 = -0.2e1 * (t258 * t357 - t278 * t358) / t272 ^ 2;
	t284 = t293 * t314 + t324 * t311;
	t310 = sin(qJ(5));
	t313 = cos(qJ(5));
	t267 = t284 * t313 + t292 * t310;
	t263 = 0.1e1 / t267;
	t264 = 0.1e1 / t267 ^ 2;
	t287 = t292 * qJD(3);
	t261 = -t283 * qJD(4) - t287 * t314;
	t288 = t293 * qJD(3);
	t253 = t267 * qJD(5) + t261 * t310 - t288 * t313;
	t266 = t284 * t310 - t292 * t313;
	t262 = t266 ^ 2;
	t257 = t262 * t264 + 0.1e1;
	t362 = t264 * t266;
	t351 = qJD(5) * t266;
	t254 = t261 * t313 + t288 * t310 - t351;
	t365 = t254 * t263 * t264;
	t367 = (t253 * t362 - t262 * t365) / t257 ^ 2;
	t366 = t250 * t283;
	t260 = t284 * qJD(4) - t287 * t311;
	t364 = t260 * t250;
	t363 = t263 * t310;
	t361 = t266 * t313;
	t360 = t268 * t283;
	t359 = t269 * t283;
	t356 = t292 * t311;
	t355 = t292 * t314;
	t352 = qJD(4) * t314;
	t279 = t283 ^ 2;
	t248 = t279 * t250 + 0.1e1;
	t350 = 0.2e1 * (-t279 * t380 + t283 * t364) / t248 ^ 2;
	t349 = -0.2e1 * t367;
	t347 = t266 * t365;
	t346 = -0.2e1 * t280 * t358;
	t345 = qJD(5) * t355 - t287;
	t334 = t264 * t361 - t363;
	t333 = -t282 * t294 + t298 * t357;
	t330 = qJD(4) * t356 + qJD(5) * t293 - t288 * t314;
	t300 = t304 * qJD(3);
	t286 = t291 * qJD(3);
	t277 = -t297 * qJD(4) + t299 * t314;
	t275 = t293 * t310 - t313 * t355;
	t274 = -t293 * t313 - t310 * t355;
	t259 = -t280 * qJD(4) + t285 * t314;
	t255 = 0.1e1 / t257;
	t246 = 0.1e1 / t248;
	t244 = t270 * t378;
	t243 = t333 * t270;
	t240 = (-t268 * t290 + t269 * t303) * t311 + t336 * t244;
	t239 = t336 * t243 - t268 * t282 + t269 * t298;
	t236 = t333 * t376 + (t298 * t346 - t259 * t294 + (t258 * t298 + t276 * t282 + t277 * t280) * t295) * t270;
	t235 = t376 * t378 + (t332 * t352 + (t303 * t346 + t286 * t294 + (t258 * t303 + t276 * t290 - t280 * t300) * t295) * t311) * t270;
	t1 = [0, 0, t235, t236, 0, 0; 0, 0, (t240 * t366 + t249 * t356) * t350 + ((-t288 * t311 - t292 * t352) * t249 + (-t364 + t379) * t240 + (t356 * t238 - (t303 * t352 - t235 * t280 - t244 * t258 - t300 * t311 + (-t244 * t297 - t290 * t311) * t242) * t359 - (-t290 * t352 - t235 * t297 - t244 * t276 + t286 * t311 + (t244 * t280 - t303 * t311) * t242) * t360) * t250) * t246, (t239 * t366 - t249 * t284) * t350 + (t239 * t379 + t261 * t249 + (-t284 * t238 - t239 * t260 - (-t236 * t280 - t243 * t258 + t277 + (-t243 * t297 - t282) * t242) * t359 - (-t236 * t297 - t243 * t276 - t259 + (t243 * t280 - t298) * t242) * t360) * t250) * t246, 0, 0; 0, 0, 0.2e1 * (-t263 * t274 + t275 * t362) * t367 + (0.2e1 * t275 * t347 - t345 * t263 * t313 + t330 * t363 + (-t345 * t266 * t310 - t275 * t253 - t274 * t254 - t330 * t361) * t264) * t255, t334 * t283 * t349 + (t334 * t260 + ((-qJD(5) * t263 - 0.2e1 * t347) * t313 + (t253 * t313 + (t254 - t351) * t310) * t264) * t283) * t255, t349 + 0.2e1 * (t253 * t264 * t255 + (-t255 * t365 - t264 * t367) * t266) * t266, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:16:22
	% EndTime: 2019-10-09 21:16:27
	% DurationCPUTime: 4.71s
	% Computational Cost: add. (21168->168), mult. (62476->326), div. (784->12), fcn. (81881->17), ass. (0->139)
	t416 = sin(pkin(12));
	t417 = sin(pkin(11));
	t383 = t417 * t416;
	t420 = cos(pkin(12));
	t421 = cos(pkin(11));
	t389 = t421 * t420;
	t423 = cos(pkin(6));
	t373 = -t423 * t389 + t383;
	t418 = sin(pkin(7));
	t419 = sin(pkin(6));
	t387 = t419 * t418;
	t422 = cos(pkin(7));
	t433 = -t373 * t422 - t421 * t387;
	t354 = sin(qJ(3));
	t385 = t417 * t420;
	t388 = t421 * t416;
	t374 = t423 * t388 + t385;
	t424 = cos(qJ(3));
	t334 = t433 * t354 + t374 * t424;
	t390 = t422 * t419;
	t343 = t373 * t418 - t421 * t390;
	t353 = sin(qJ(4));
	t356 = cos(qJ(4));
	t320 = -t334 * t353 + t343 * t356;
	t428 = -t374 * t354 + t433 * t424;
	t328 = t428 * qJD(3);
	t298 = t320 * qJD(4) + t328 * t356;
	t321 = t334 * t356 + t343 * t353;
	t355 = cos(qJ(5));
	t352 = sin(qJ(5));
	t429 = t428 * t352;
	t308 = t321 * t355 - t429;
	t363 = t334 * qJD(3);
	t280 = t308 * qJD(5) + t298 * t352 - t355 * t363;
	t364 = t428 * t355;
	t306 = t321 * t352 + t364;
	t301 = t306 ^ 2;
	t386 = t419 * t416;
	t426 = t420 * t390 + t423 * t418;
	t342 = t426 * t354 + t424 * t386;
	t347 = -t420 * t387 + t423 * t422;
	t338 = t342 * t356 + t347 * t353;
	t369 = -t354 * t386 + t426 * t424;
	t324 = t338 * t352 + t355 * t369;
	t318 = 0.1e1 / t324 ^ 2;
	t292 = t301 * t318 + 0.1e1;
	t289 = 0.1e1 / t292;
	t337 = -t342 * t353 + t347 * t356;
	t339 = t369 * qJD(3);
	t315 = t337 * qJD(4) + t339 * t356;
	t325 = t338 * t355 - t352 * t369;
	t340 = t342 * qJD(3);
	t294 = t325 * qJD(5) + t315 * t352 - t340 * t355;
	t317 = 0.1e1 / t324;
	t407 = t306 * t318;
	t267 = (-t280 * t317 + t294 * t407) * t289;
	t293 = atan2(-t306, t324);
	t285 = sin(t293);
	t286 = cos(t293);
	t382 = -t285 * t324 - t286 * t306;
	t263 = t382 * t267 - t285 * t280 + t286 * t294;
	t279 = -t285 * t306 + t286 * t324;
	t276 = 0.1e1 / t279;
	t277 = 0.1e1 / t279 ^ 2;
	t432 = t263 * t276 * t277;
	t348 = -t423 * t383 + t389;
	t375 = t423 * t385 + t388;
	t384 = t417 * t419;
	t427 = t375 * t422 - t418 * t384;
	t336 = t348 * t424 - t427 * t354;
	t344 = t375 * t418 + t422 * t384;
	t323 = t336 * t356 + t344 * t353;
	t335 = t348 * t354 + t427 * t424;
	t309 = t323 * t352 - t335 * t355;
	t431 = -0.2e1 * t309;
	t392 = 0.2e1 * t309 * t432;
	t379 = -t317 * t320 + t337 * t407;
	t430 = t352 * t379;
	t409 = t294 * t317 * t318;
	t425 = -0.2e1 * (t280 * t407 - t301 * t409) / t292 ^ 2;
	t310 = t323 * t355 + t335 * t352;
	t303 = 0.1e1 / t310;
	t304 = 0.1e1 / t310 ^ 2;
	t322 = -t336 * t353 + t344 * t356;
	t316 = t322 ^ 2;
	t406 = t316 * t304;
	t291 = 0.1e1 + t406;
	t329 = t335 * qJD(3);
	t299 = -t323 * qJD(4) + t329 * t353;
	t300 = t322 * qJD(4) - t329 * t356;
	t330 = t336 * qJD(3);
	t283 = -t309 * qJD(5) + t300 * t355 + t330 * t352;
	t412 = t283 * t303 * t304;
	t395 = t316 * t412;
	t408 = t304 * t322;
	t415 = (t299 * t408 - t395) / t291 ^ 2;
	t414 = t277 * t309;
	t282 = t310 * qJD(5) + t300 * t352 - t330 * t355;
	t413 = t282 * t277;
	t411 = t285 * t309;
	t410 = t286 * t309;
	t405 = t322 * t352;
	t404 = t352 * t356;
	t403 = t355 * t356;
	t402 = qJD(4) * t353;
	t401 = qJD(5) * t352;
	t400 = qJD(5) * t355;
	t399 = t356 * qJD(5);
	t302 = t309 ^ 2;
	t275 = t302 * t277 + 0.1e1;
	t398 = 0.2e1 * (-t302 * t432 + t309 * t413) / t275 ^ 2;
	t397 = 0.2e1 * t415;
	t394 = t322 * t412;
	t393 = -0.2e1 * t306 * t409;
	t381 = -t308 * t317 + t325 * t407;
	t365 = t356 * t428;
	t311 = -t334 * t355 + t352 * t365;
	t326 = -t342 * t355 + t369 * t404;
	t380 = -t311 * t317 + t326 * t407;
	t314 = -t338 * qJD(4) - t339 * t353;
	t313 = -t335 * t403 + t336 * t352;
	t312 = -t335 * t404 - t336 * t355;
	t297 = -t321 * qJD(4) - t328 * t353;
	t296 = (t369 * t399 - t339) * t355 + (qJD(5) * t342 - t340 * t356 - t369 * t402) * t352;
	t295 = -t324 * qJD(5) + t315 * t355 + t340 * t352;
	t287 = 0.1e1 / t291;
	t284 = -t328 * t355 + t334 * t401 - t404 * t363 + t365 * t400 - t402 * t429;
	t281 = -qJD(5) * t364 + t298 * t355 - t321 * t401 + t352 * t363;
	t273 = 0.1e1 / t275;
	t272 = t289 * t430;
	t271 = t380 * t289;
	t269 = t381 * t289;
	t266 = (-t285 * t320 + t286 * t337) * t352 + t382 * t272;
	t265 = t382 * t271 - t285 * t311 + t286 * t326;
	t264 = t382 * t269 - t285 * t308 + t286 * t325;
	t262 = t380 * t425 + (t326 * t393 - t284 * t317 + (t280 * t326 + t294 * t311 + t296 * t306) * t318) * t289;
	t260 = t381 * t425 + (t325 * t393 - t281 * t317 + (t280 * t325 + t294 * t308 + t295 * t306) * t318) * t289;
	t259 = t425 * t430 + (t379 * t400 + (t337 * t393 - t297 * t317 + (t280 * t337 + t294 * t320 + t306 * t314) * t318) * t352) * t289;
	t1 = [0, 0, t262, t259, t260, 0; 0, 0, (t265 * t414 - t276 * t312) * t398 + (t265 * t392 + (-t312 * t263 - t265 * t282 - (-t262 * t306 - t271 * t280 + t296 + (-t271 * t324 - t311) * t267) * t410 - (-t262 * t324 - t271 * t294 - t284 + (t271 * t306 - t326) * t267) * t411) * t277 + ((-t335 * t399 + t329) * t355 + (qJD(5) * t336 - t330 * t356 + t335 * t402) * t352) * t276) * t273, (t266 * t414 - t276 * t405) * t398 + ((t299 * t352 + t322 * t400) * t276 + (-t413 + t392) * t266 + (-t405 * t263 - (t337 * t400 - t259 * t306 - t272 * t280 + t314 * t352 + (-t272 * t324 - t320 * t352) * t267) * t410 - (-t320 * t400 - t259 * t324 - t272 * t294 - t297 * t352 + (t272 * t306 - t337 * t352) * t267) * t411) * t277) * t273, (t264 * t414 - t276 * t310) * t398 + (t264 * t392 + t283 * t276 + (-t310 * t263 - t264 * t282 - (-t260 * t306 - t269 * t280 + t295 + (-t269 * t324 - t308) * t267) * t410 - (-t260 * t324 - t269 * t294 - t281 + (t269 * t306 - t325) * t267) * t411) * t277) * t273, 0; 0, 0, (-t303 * t335 * t353 + t313 * t408) * t397 + (0.2e1 * t313 * t394 + (t335 * qJD(4) * t356 + t330 * t353) * t303 + (-(-t329 * t352 - t330 * t403 + t336 * t400) * t322 - t313 * t299 + (-t353 * t283 - (t352 * t399 + t355 * t402) * t322) * t335) * t304) * t287, (t303 * t323 + t355 * t406) * t397 + (0.2e1 * t355 * t395 - t300 * t303 + (-0.2e1 * t299 * t322 * t355 + t283 * t323 + t316 * t401) * t304) * t287, t408 * t415 * t431 + (t394 * t431 + (t282 * t322 + t299 * t309) * t304) * t287, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end