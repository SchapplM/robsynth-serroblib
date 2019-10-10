% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR14_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR14_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiR_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t50 = sin(qJ(2));
	t51 = sin(qJ(1));
	t57 = t51 * t50;
	t52 = cos(qJ(2));
	t56 = t51 * t52;
	t53 = cos(qJ(1));
	t55 = t53 * t50;
	t54 = t53 * t52;
	t49 = cos(pkin(6));
	t48 = sin(pkin(6));
	t47 = -t49 * t57 + t54;
	t46 = -t49 * t56 - t55;
	t45 = -t49 * t55 - t56;
	t44 = -t49 * t54 + t57;
	t1 = [t45, t46, 0, 0, 0, 0; t47, -t44, 0, 0, 0, 0; 0, t48 * t52, 0, 0, 0, 0; t44, -t47, 0, 0, 0, 0; t46, t45, 0, 0, 0, 0; 0, -t48 * t50, 0, 0, 0, 0; t53 * t48, 0, 0, 0, 0, 0; t51 * t48, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (38->19), mult. (116->49), div. (0->0), fcn. (167->10), ass. (0->24)
	t86 = sin(pkin(14));
	t90 = cos(pkin(7));
	t104 = t86 * t90;
	t88 = sin(pkin(6));
	t92 = sin(qJ(1));
	t103 = t88 * t92;
	t94 = cos(qJ(1));
	t102 = t88 * t94;
	t89 = cos(pkin(14));
	t101 = t89 * t90;
	t91 = sin(qJ(2));
	t100 = t90 * t91;
	t99 = cos(pkin(6));
	t98 = t92 * t99;
	t97 = t94 * t99;
	t93 = cos(qJ(2));
	t80 = t92 * t91 - t93 * t97;
	t87 = sin(pkin(7));
	t96 = t87 * t102 + t80 * t90;
	t82 = -t94 * t91 - t93 * t98;
	t95 = t87 * t103 + t82 * t90;
	t83 = -t91 * t98 + t94 * t93;
	t81 = -t91 * t97 - t92 * t93;
	t1 = [t81 * t89 + t96 * t86, -t83 * t104 + t82 * t89, 0, 0, 0, 0; t83 * t89 + t95 * t86, t81 * t104 - t80 * t89, 0, 0, 0, 0; 0, (-t86 * t100 + t89 * t93) * t88, 0, 0, 0, 0; -t81 * t86 + t96 * t89, -t83 * t101 - t82 * t86, 0, 0, 0, 0; -t83 * t86 + t95 * t89, t81 * t101 + t80 * t86, 0, 0, 0, 0; 0, (-t89 * t100 - t86 * t93) * t88, 0, 0, 0, 0; t90 * t102 - t80 * t87, t83 * t87, 0, 0, 0, 0; t90 * t103 - t82 * t87, -t81 * t87, 0, 0, 0, 0; 0, t88 * t91 * t87, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:59
	% EndTime: 2019-10-10 11:10:59
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (180->48), mult. (545->107), div. (0->0), fcn. (746->14), ass. (0->57)
	t199 = sin(qJ(2));
	t200 = sin(qJ(1));
	t202 = cos(qJ(2));
	t203 = cos(qJ(1));
	t228 = cos(pkin(6));
	t214 = t203 * t228;
	t186 = t199 * t214 + t200 * t202;
	t191 = sin(pkin(14));
	t195 = cos(pkin(14));
	t185 = t200 * t199 - t202 * t214;
	t193 = sin(pkin(7));
	t197 = cos(pkin(7));
	t194 = sin(pkin(6));
	t220 = t194 * t203;
	t206 = t185 * t197 + t193 * t220;
	t170 = -t186 * t195 + t206 * t191;
	t198 = sin(qJ(4));
	t201 = cos(qJ(4));
	t169 = t186 * t191 + t206 * t195;
	t179 = -t185 * t193 + t197 * t220;
	t192 = sin(pkin(8));
	t196 = cos(pkin(8));
	t211 = t169 * t196 + t179 * t192;
	t234 = t170 * t201 + t211 * t198;
	t233 = -t170 * t198 + t211 * t201;
	t225 = t191 * t197;
	t224 = t192 * t193;
	t223 = t193 * t196;
	t222 = t194 * t199;
	t221 = t194 * t200;
	t219 = t195 * t197;
	t218 = t197 * t199;
	t217 = t197 * t202;
	t216 = t193 * t222;
	t215 = t200 * t228;
	t213 = t228 * t193;
	t188 = -t199 * t215 + t203 * t202;
	t187 = -t203 * t199 - t202 * t215;
	t205 = t187 * t197 + t193 * t221;
	t171 = -t188 * t191 + t205 * t195;
	t181 = -t187 * t193 + t197 * t221;
	t210 = t171 * t196 + t181 * t192;
	t209 = (t195 * t213 + (-t191 * t199 + t195 * t217) * t194) * t196 + (-t194 * t202 * t193 + t228 * t197) * t192;
	t173 = t185 * t191 - t186 * t219;
	t208 = t173 * t196 + t186 * t224;
	t175 = -t187 * t191 - t188 * t219;
	t207 = t175 * t196 + t188 * t224;
	t182 = (-t191 * t202 - t195 * t218) * t194;
	t204 = t182 * t196 + t192 * t216;
	t183 = (-t191 * t218 + t195 * t202) * t194;
	t178 = t195 * t222 + (t194 * t217 + t213) * t191;
	t176 = t187 * t195 - t188 * t225;
	t174 = -t185 * t195 - t186 * t225;
	t172 = t188 * t195 + t205 * t191;
	t166 = t172 * t201 + t210 * t198;
	t165 = -t172 * t198 + t210 * t201;
	t1 = [t234, t176 * t201 + t207 * t198, 0, t165, 0, 0; t166, t174 * t201 + t208 * t198, 0, -t233, 0, 0; 0, t183 * t201 + t204 * t198, 0, -t178 * t198 + t209 * t201, 0, 0; t233, -t176 * t198 + t207 * t201, 0, -t166, 0, 0; t165, -t174 * t198 + t208 * t201, 0, t234, 0, 0; 0, -t183 * t198 + t204 * t201, 0, -t178 * t201 - t209 * t198, 0, 0; -t169 * t192 + t179 * t196, -t175 * t192 + t188 * t223, 0, 0, 0, 0; -t171 * t192 + t181 * t196, -t173 * t192 + t186 * t223, 0, 0, 0, 0; 0, -t182 * t192 + t196 * t216, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:11:00
	% EndTime: 2019-10-10 11:11:01
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (400->65), mult. (1183->139), div. (0->0), fcn. (1610->16), ass. (0->76)
	t298 = sin(qJ(2));
	t299 = sin(qJ(1));
	t302 = cos(qJ(2));
	t303 = cos(qJ(1));
	t329 = cos(pkin(6));
	t313 = t303 * t329;
	t284 = t298 * t313 + t299 * t302;
	t289 = sin(pkin(14));
	t293 = cos(pkin(14));
	t283 = t299 * t298 - t302 * t313;
	t291 = sin(pkin(7));
	t295 = cos(pkin(7));
	t292 = sin(pkin(6));
	t319 = t292 * t303;
	t306 = t283 * t295 + t291 * t319;
	t266 = -t284 * t293 + t306 * t289;
	t297 = sin(qJ(4));
	t301 = cos(qJ(4));
	t265 = t284 * t289 + t306 * t293;
	t277 = -t283 * t291 + t295 * t319;
	t290 = sin(pkin(8));
	t294 = cos(pkin(8));
	t310 = t265 * t294 + t277 * t290;
	t248 = t266 * t301 + t310 * t297;
	t255 = t265 * t290 - t277 * t294;
	t296 = sin(qJ(5));
	t300 = cos(qJ(5));
	t338 = t248 * t296 + t255 * t300;
	t337 = t248 * t300 - t255 * t296;
	t246 = t266 * t297 - t310 * t301;
	t314 = t299 * t329;
	t286 = -t298 * t314 + t303 * t302;
	t285 = -t303 * t298 - t302 * t314;
	t320 = t292 * t299;
	t305 = t285 * t295 + t291 * t320;
	t267 = -t286 * t289 + t305 * t293;
	t279 = -t285 * t291 + t295 * t320;
	t330 = t267 * t294 + t279 * t290;
	t324 = t289 * t295;
	t323 = t290 * t291;
	t322 = t291 * t294;
	t321 = t292 * t298;
	t318 = t293 * t295;
	t317 = t295 * t298;
	t316 = t295 * t302;
	t315 = t291 * t321;
	t312 = t329 * t291;
	t275 = t293 * t312 + (-t289 * t298 + t293 * t316) * t292;
	t282 = -t292 * t302 * t291 + t329 * t295;
	t309 = t275 * t294 + t282 * t290;
	t269 = t283 * t289 - t284 * t318;
	t308 = t269 * t294 + t284 * t323;
	t271 = -t285 * t289 - t286 * t318;
	t307 = t271 * t294 + t286 * t323;
	t280 = (-t289 * t302 - t293 * t317) * t292;
	t304 = t280 * t294 + t290 * t315;
	t281 = (-t289 * t317 + t293 * t302) * t292;
	t276 = t293 * t321 + (t292 * t316 + t312) * t289;
	t273 = -t280 * t290 + t294 * t315;
	t272 = t285 * t293 - t286 * t324;
	t270 = -t283 * t293 - t284 * t324;
	t268 = t286 * t293 + t305 * t289;
	t262 = -t275 * t290 + t282 * t294;
	t260 = -t271 * t290 + t286 * t322;
	t259 = -t269 * t290 + t284 * t322;
	t258 = t281 * t301 + t304 * t297;
	t257 = -t267 * t290 + t279 * t294;
	t254 = t276 * t301 + t309 * t297;
	t253 = -t276 * t297 + t309 * t301;
	t252 = t272 * t301 + t307 * t297;
	t251 = t270 * t301 + t308 * t297;
	t250 = t268 * t301 + t330 * t297;
	t249 = t268 * t297 - t330 * t301;
	t245 = t250 * t300 + t257 * t296;
	t244 = -t250 * t296 + t257 * t300;
	t1 = [t337, t252 * t300 + t260 * t296, 0, -t249 * t300, t244, 0; t245, t251 * t300 + t259 * t296, 0, t246 * t300, t338, 0; 0, t258 * t300 + t273 * t296, 0, t253 * t300, -t254 * t296 + t262 * t300, 0; -t338, -t252 * t296 + t260 * t300, 0, t249 * t296, -t245, 0; t244, -t251 * t296 + t259 * t300, 0, -t246 * t296, t337, 0; 0, -t258 * t296 + t273 * t300, 0, -t253 * t296, -t254 * t300 - t262 * t296, 0; t246, t272 * t297 - t307 * t301, 0, t250, 0, 0; t249, t270 * t297 - t308 * t301, 0, -t248, 0, 0; 0, t281 * t297 - t304 * t301, 0, t254, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:11:01
	% EndTime: 2019-10-10 11:11:03
	% DurationCPUTime: 1.05s
	% Computational Cost: add. (825->85), mult. (2425->180), div. (0->0), fcn. (3286->18), ass. (0->92)
	t368 = sin(qJ(2));
	t369 = sin(qJ(1));
	t373 = cos(qJ(2));
	t374 = cos(qJ(1));
	t412 = cos(pkin(6));
	t394 = t374 * t412;
	t353 = t368 * t394 + t369 * t373;
	t358 = sin(pkin(14));
	t362 = cos(pkin(14));
	t352 = t369 * t368 - t373 * t394;
	t360 = sin(pkin(7));
	t364 = cos(pkin(7));
	t361 = sin(pkin(6));
	t402 = t361 * t374;
	t389 = t352 * t364 + t360 * t402;
	t337 = -t353 * t362 + t389 * t358;
	t367 = sin(qJ(4));
	t372 = cos(qJ(4));
	t348 = -t352 * t360 + t364 * t402;
	t359 = sin(pkin(8));
	t363 = cos(pkin(8));
	t413 = t353 * t358 + t389 * t362;
	t417 = t348 * t359 + t363 * t413;
	t315 = t337 * t372 + t367 * t417;
	t326 = -t348 * t363 + t359 * t413;
	t366 = sin(qJ(5));
	t371 = cos(qJ(5));
	t305 = t315 * t371 - t326 * t366;
	t314 = t337 * t367 - t372 * t417;
	t365 = sin(qJ(6));
	t370 = cos(qJ(6));
	t423 = t305 * t365 - t314 * t370;
	t422 = t305 * t370 + t314 * t365;
	t303 = t315 * t366 + t326 * t371;
	t393 = t412 * t360;
	t399 = t364 * t373;
	t377 = t362 * t393 + (-t358 * t368 + t362 * t399) * t361;
	t385 = -t361 * t373 * t360 + t412 * t364;
	t416 = t385 * t359 + t377 * t363;
	t395 = t369 * t412;
	t355 = -t368 * t395 + t374 * t373;
	t354 = -t374 * t368 - t373 * t395;
	t403 = t361 * t369;
	t387 = t354 * t364 + t360 * t403;
	t380 = t355 * t358 - t387 * t362;
	t388 = -t354 * t360 + t364 * t403;
	t415 = -t388 * t359 + t380 * t363;
	t407 = t358 * t364;
	t406 = t360 * t359;
	t405 = t360 * t363;
	t404 = t361 * t368;
	t401 = t362 * t364;
	t400 = t364 * t368;
	t398 = t365 * t371;
	t397 = t370 * t371;
	t396 = t360 * t404;
	t339 = t352 * t358 - t353 * t401;
	t391 = t339 * t363 + t353 * t406;
	t341 = -t354 * t358 - t355 * t401;
	t390 = t341 * t363 + t355 * t406;
	t350 = (-t358 * t373 - t362 * t400) * t361;
	t386 = t350 * t363 + t359 * t396;
	t375 = t380 * t359 + t388 * t363;
	t351 = (-t358 * t400 + t362 * t373) * t361;
	t347 = t362 * t404 + (t361 * t399 + t393) * t358;
	t343 = -t350 * t359 + t363 * t396;
	t342 = t354 * t362 - t355 * t407;
	t340 = -t352 * t362 - t353 * t407;
	t338 = t355 * t362 + t387 * t358;
	t334 = -t377 * t359 + t385 * t363;
	t331 = -t341 * t359 + t355 * t405;
	t330 = -t339 * t359 + t353 * t405;
	t329 = t351 * t372 + t367 * t386;
	t328 = t351 * t367 - t372 * t386;
	t324 = t347 * t372 + t416 * t367;
	t323 = t347 * t367 - t416 * t372;
	t322 = t329 * t371 + t343 * t366;
	t321 = t342 * t372 + t367 * t390;
	t320 = t342 * t367 - t372 * t390;
	t319 = t340 * t372 + t367 * t391;
	t318 = t340 * t367 - t372 * t391;
	t317 = t338 * t372 - t415 * t367;
	t316 = t338 * t367 + t415 * t372;
	t311 = t324 * t371 + t334 * t366;
	t310 = -t324 * t366 + t334 * t371;
	t309 = t321 * t371 + t331 * t366;
	t308 = t319 * t371 + t330 * t366;
	t307 = t317 * t371 + t366 * t375;
	t306 = t317 * t366 - t371 * t375;
	t302 = t307 * t370 + t316 * t365;
	t301 = -t307 * t365 + t316 * t370;
	t1 = [t422, t309 * t370 + t320 * t365, 0, -t316 * t397 + t317 * t365, -t306 * t370, t301; t302, t308 * t370 + t318 * t365, 0, t314 * t397 - t315 * t365, t303 * t370, t423; 0, t322 * t370 + t328 * t365, 0, -t323 * t397 + t324 * t365, t310 * t370, -t311 * t365 + t323 * t370; -t423, -t309 * t365 + t320 * t370, 0, t316 * t398 + t317 * t370, t306 * t365, -t302; t301, -t308 * t365 + t318 * t370, 0, -t314 * t398 - t315 * t370, -t303 * t365, t422; 0, -t322 * t365 + t328 * t370, 0, t323 * t398 + t324 * t370, -t310 * t365, -t311 * t370 - t323 * t365; t303, t321 * t366 - t331 * t371, 0, -t316 * t366, t307, 0; t306, t319 * t366 - t330 * t371, 0, t314 * t366, -t305, 0; 0, t329 * t366 - t343 * t371, 0, -t323 * t366, t311, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end