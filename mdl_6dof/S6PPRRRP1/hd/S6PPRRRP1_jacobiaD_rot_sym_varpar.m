% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRRP1
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
%   Wie in S6PPRRRP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:14
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRRP1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRP1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:14:28
	% EndTime: 2019-10-09 21:14:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:14:28
	% EndTime: 2019-10-09 21:14:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:14:28
	% EndTime: 2019-10-09 21:14:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:14:28
	% EndTime: 2019-10-09 21:14:28
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
	% StartTime: 2019-10-09 21:14:29
	% EndTime: 2019-10-09 21:14:29
	% DurationCPUTime: 0.96s
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
	% StartTime: 2019-10-09 21:14:29
	% EndTime: 2019-10-09 21:14:31
	% DurationCPUTime: 2.51s
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
	% StartTime: 2019-10-09 21:14:29
	% EndTime: 2019-10-09 21:14:32
	% DurationCPUTime: 2.51s
	% Computational Cost: add. (9293->119), mult. (27507->244), div. (559->12), fcn. (36133->17), ass. (0->119)
	t379 = sin(pkin(12));
	t380 = sin(pkin(11));
	t348 = t380 * t379;
	t383 = cos(pkin(12));
	t384 = cos(pkin(11));
	t354 = t384 * t383;
	t386 = cos(pkin(6));
	t319 = -t386 * t348 + t354;
	t323 = sin(qJ(3));
	t326 = cos(qJ(3));
	t350 = t380 * t383;
	t353 = t384 * t379;
	t340 = t386 * t350 + t353;
	t382 = sin(pkin(6));
	t349 = t380 * t382;
	t381 = sin(pkin(7));
	t385 = cos(pkin(7));
	t388 = t340 * t385 - t381 * t349;
	t303 = t319 * t323 + t388 * t326;
	t318 = t386 * t353 + t350;
	t339 = -t386 * t354 + t348;
	t352 = t382 * t381;
	t334 = -t339 * t385 - t384 * t352;
	t302 = t318 * t326 + t334 * t323;
	t322 = sin(qJ(4));
	t325 = cos(qJ(4));
	t355 = t385 * t382;
	t333 = t339 * t381 - t384 * t355;
	t293 = t302 * t325 + t333 * t322;
	t301 = -t318 * t323 + t334 * t326;
	t296 = t301 * qJD(3);
	t269 = t293 * qJD(4) + t296 * t322;
	t291 = t302 * t322 - t333 * t325;
	t289 = t291 ^ 2;
	t338 = t383 * t355 + t386 * t381;
	t351 = t382 * t379;
	t315 = t338 * t323 + t326 * t351;
	t317 = -t383 * t352 + t386 * t385;
	t308 = t315 * t322 - t317 * t325;
	t306 = 0.1e1 / t308 ^ 2;
	t283 = t289 * t306 + 0.1e1;
	t281 = 0.1e1 / t283;
	t309 = t315 * t325 + t317 * t322;
	t314 = -t323 * t351 + t338 * t326;
	t310 = t314 * qJD(3);
	t287 = t309 * qJD(4) + t310 * t322;
	t305 = 0.1e1 / t308;
	t368 = t291 * t306;
	t253 = (-t269 * t305 + t287 * t368) * t281;
	t284 = atan2(-t291, t308);
	t279 = sin(t284);
	t280 = cos(t284);
	t347 = -t279 * t308 - t280 * t291;
	t249 = t347 * t253 - t279 * t269 + t280 * t287;
	t263 = -t279 * t291 + t280 * t308;
	t260 = 0.1e1 / t263;
	t261 = 0.1e1 / t263 ^ 2;
	t391 = t249 * t260 * t261;
	t304 = t319 * t326 - t323 * t388;
	t335 = t340 * t381 + t385 * t349;
	t294 = t304 * t322 - t335 * t325;
	t390 = 0.2e1 * t294 * t391;
	t343 = -t301 * t305 + t314 * t368;
	t389 = t322 * t343;
	t369 = t287 * t305 * t306;
	t387 = -0.2e1 * (t269 * t368 - t289 * t369) / t283 ^ 2;
	t295 = t304 * t325 + t335 * t322;
	t321 = sin(qJ(5));
	t324 = cos(qJ(5));
	t278 = t295 * t324 + t303 * t321;
	t274 = 0.1e1 / t278;
	t275 = 0.1e1 / t278 ^ 2;
	t298 = t303 * qJD(3);
	t272 = -t294 * qJD(4) - t298 * t325;
	t299 = t304 * qJD(3);
	t264 = t278 * qJD(5) + t272 * t321 - t299 * t324;
	t277 = t295 * t321 - t303 * t324;
	t273 = t277 ^ 2;
	t268 = t273 * t275 + 0.1e1;
	t373 = t275 * t277;
	t362 = qJD(5) * t277;
	t265 = t272 * t324 + t299 * t321 - t362;
	t376 = t265 * t274 * t275;
	t378 = (t264 * t373 - t273 * t376) / t268 ^ 2;
	t377 = t261 * t294;
	t271 = t295 * qJD(4) - t298 * t322;
	t375 = t271 * t261;
	t374 = t274 * t321;
	t372 = t277 * t324;
	t371 = t279 * t294;
	t370 = t280 * t294;
	t367 = t303 * t322;
	t366 = t303 * t325;
	t363 = qJD(4) * t325;
	t290 = t294 ^ 2;
	t259 = t290 * t261 + 0.1e1;
	t361 = 0.2e1 * (-t290 * t391 + t294 * t375) / t259 ^ 2;
	t360 = -0.2e1 * t378;
	t358 = t277 * t376;
	t357 = -0.2e1 * t291 * t369;
	t356 = qJD(5) * t366 - t298;
	t345 = t275 * t372 - t374;
	t344 = -t293 * t305 + t309 * t368;
	t341 = qJD(4) * t367 + qJD(5) * t304 - t299 * t325;
	t311 = t315 * qJD(3);
	t297 = t302 * qJD(3);
	t288 = -t308 * qJD(4) + t310 * t325;
	t286 = t304 * t321 - t324 * t366;
	t285 = -t304 * t324 - t321 * t366;
	t270 = -t291 * qJD(4) + t296 * t325;
	t266 = 0.1e1 / t268;
	t257 = 0.1e1 / t259;
	t255 = t281 * t389;
	t254 = t344 * t281;
	t251 = (-t279 * t301 + t280 * t314) * t322 + t347 * t255;
	t250 = t347 * t254 - t279 * t293 + t280 * t309;
	t247 = t344 * t387 + (t309 * t357 - t270 * t305 + (t269 * t309 + t287 * t293 + t288 * t291) * t306) * t281;
	t246 = t387 * t389 + (t343 * t363 + (t314 * t357 + t297 * t305 + (t269 * t314 + t287 * t301 - t291 * t311) * t306) * t322) * t281;
	t1 = [0, 0, t246, t247, 0, 0; 0, 0, (t251 * t377 + t260 * t367) * t361 + ((-t299 * t322 - t303 * t363) * t260 + (-t375 + t390) * t251 + (t367 * t249 - (t314 * t363 - t246 * t291 - t255 * t269 - t311 * t322 + (-t255 * t308 - t301 * t322) * t253) * t370 - (-t301 * t363 - t246 * t308 - t255 * t287 + t297 * t322 + (t255 * t291 - t314 * t322) * t253) * t371) * t261) * t257, (t250 * t377 - t260 * t295) * t361 + (t250 * t390 + t272 * t260 + (-t295 * t249 - t250 * t271 - (-t247 * t291 - t254 * t269 + t288 + (-t254 * t308 - t293) * t253) * t370 - (-t247 * t308 - t254 * t287 - t270 + (t254 * t291 - t309) * t253) * t371) * t261) * t257, 0, 0; 0, 0, 0.2e1 * (-t274 * t285 + t286 * t373) * t378 + (0.2e1 * t286 * t358 - t356 * t274 * t324 + t341 * t374 + (-t356 * t277 * t321 - t286 * t264 - t285 * t265 - t341 * t372) * t275) * t266, t345 * t294 * t360 + (t345 * t271 + ((-qJD(5) * t274 - 0.2e1 * t358) * t324 + (t264 * t324 + (t265 - t362) * t321) * t275) * t294) * t266, t360 + 0.2e1 * (t264 * t275 * t266 + (-t266 * t376 - t275 * t378) * t277) * t277, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end