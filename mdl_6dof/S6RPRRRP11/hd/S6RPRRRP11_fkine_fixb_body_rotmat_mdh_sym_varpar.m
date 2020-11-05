% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRP11 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:54
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRRP11_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:54:11
	% EndTime: 2020-11-04 21:54:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:54:11
	% EndTime: 2020-11-04 21:54:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t172 = cos(qJ(1));
	t171 = sin(qJ(1));
	t1 = [t172, -t171, 0, 0; t171, t172, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:54:11
	% EndTime: 2020-11-04 21:54:11
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t173 = sin(pkin(12));
	t177 = sin(qJ(1));
	t184 = t177 * t173;
	t174 = sin(pkin(6));
	t183 = t177 * t174;
	t175 = cos(pkin(12));
	t182 = t177 * t175;
	t178 = cos(qJ(1));
	t181 = t178 * t173;
	t180 = t178 * t174;
	t179 = t178 * t175;
	t176 = cos(pkin(6));
	t1 = [-t176 * t184 + t179, -t176 * t182 - t181, t183, t178 * pkin(1) + qJ(2) * t183 + 0; t176 * t181 + t182, t176 * t179 - t184, -t180, t177 * pkin(1) - qJ(2) * t180 + 0; t174 * t173, t174 * t175, t176, t176 * qJ(2) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:54:11
	% EndTime: 2020-11-04 21:54:11
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->32), mult. (101->58), div. (0->0), fcn. (135->10), ass. (0->31)
	t193 = sin(pkin(7));
	t216 = pkin(9) * t193;
	t192 = sin(pkin(12));
	t198 = sin(qJ(3));
	t215 = t192 * t198;
	t200 = cos(qJ(3));
	t214 = t192 * t200;
	t194 = sin(pkin(6));
	t213 = t194 * t193;
	t196 = cos(pkin(7));
	t212 = t194 * t196;
	t197 = cos(pkin(6));
	t211 = t197 * t196;
	t210 = t197 * t198;
	t209 = t197 * t200;
	t195 = cos(pkin(12));
	t208 = t198 * t195;
	t199 = sin(qJ(1));
	t207 = t199 * t192;
	t206 = t200 * t195;
	t201 = cos(qJ(1));
	t205 = t201 * t192;
	t204 = t201 * t195;
	t188 = -t192 * pkin(2) + t195 * t216;
	t190 = t196 * pkin(9) + qJ(2);
	t203 = t188 * t197 + t194 * t190;
	t185 = t195 * t211 - t213;
	t202 = t185 * t198 + t192 * t209;
	t187 = t195 * pkin(2) + t192 * t216 + pkin(1);
	t186 = t196 * t215 - t206;
	t1 = [-t201 * t186 - t202 * t199, (-t185 * t199 - t196 * t205) * t200 + t198 * (t197 * t207 - t204), (t199 * t197 * t195 + t205) * t193 + t199 * t212, t187 * t201 + t203 * t199 + 0; -t199 * t186 + t202 * t201, (t185 * t200 - t192 * t210) * t201 - t199 * (t196 * t214 + t208), -(t197 * t204 - t207) * t193 - t201 * t212, t187 * t199 - t203 * t201 + 0; t193 * t210 + (t196 * t208 + t214) * t194, t193 * t209 + (t196 * t206 - t215) * t194, -t195 * t213 + t211, -t188 * t194 + t190 * t197 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:54:11
	% EndTime: 2020-11-04 21:54:11
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (91->53), mult. (213->97), div. (0->0), fcn. (268->12), ass. (0->38)
	t236 = cos(pkin(12));
	t237 = cos(pkin(7));
	t238 = cos(pkin(6));
	t252 = t238 * t237;
	t234 = sin(pkin(7));
	t235 = sin(pkin(6));
	t255 = t235 * t234;
	t224 = t236 * t252 - t255;
	t233 = sin(pkin(12));
	t240 = sin(qJ(3));
	t243 = cos(qJ(3));
	t251 = t238 * t243;
	t218 = t224 * t240 + t233 * t251;
	t253 = t238 * t234;
	t222 = t235 * t237 + t236 * t253;
	t239 = sin(qJ(4));
	t242 = cos(qJ(4));
	t262 = t218 * t239 + t242 * t222;
	t226 = t236 * t234 * pkin(9) - t233 * pkin(2);
	t254 = t236 * t237;
	t227 = -t233 * pkin(3) + pkin(10) * t254;
	t228 = pkin(3) * t254 + t233 * pkin(10);
	t230 = t237 * pkin(9) + qJ(2);
	t261 = (pkin(3) * t255 - t228 * t238) * t240 - (pkin(10) * t255 - t227 * t238) * t243 + t226 * t238 + t235 * t230;
	t260 = t233 * t234;
	t259 = t233 * t237;
	t258 = t233 * t240;
	t257 = t233 * t243;
	t244 = cos(qJ(1));
	t256 = t233 * t244;
	t249 = t243 * t236;
	t245 = (t235 * t254 + t253) * t240 + t235 * t257;
	t241 = sin(qJ(1));
	t225 = -t237 * t258 + t249;
	t221 = t236 * t255 - t252;
	t220 = t225 * t239 - t242 * t260;
	t217 = (t236 * pkin(3) + pkin(10) * t259) * t243 + (-pkin(3) * t259 + t236 * pkin(10)) * t240 + pkin(9) * t260 + t236 * pkin(2) + pkin(1);
	t1 = [(-t218 * t241 + t244 * t225) * t242 + (t222 * t241 + t234 * t256) * t239, -t220 * t244 + t262 * t241, (t224 * t241 + t237 * t256) * t243 - t240 * (t241 * t238 * t233 - t244 * t236), t217 * t244 + t261 * t241 + 0; (t218 * t242 - t239 * t222) * t244 + (t225 * t242 + t239 * t260) * t241, -t220 * t241 - t262 * t244, (-t224 * t243 + t238 * t258) * t244 + t241 * (t240 * t236 + t237 * t257), t217 * t241 - t261 * t244 + 0; -t239 * t221 + t245 * t242, -t242 * t221 - t245 * t239, -t234 * t251 + (-t237 * t249 + t258) * t235, (-pkin(10) * t253 - t227 * t235) * t243 + (pkin(3) * t253 + t228 * t235) * t240 - t226 * t235 + t230 * t238 + 0 + pkin(8); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:54:11
	% EndTime: 2020-11-04 21:54:11
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (167->66), mult. (417->110), div. (0->0), fcn. (518->14), ass. (0->53)
	t288 = cos(pkin(12));
	t289 = cos(pkin(7));
	t290 = cos(pkin(6));
	t313 = t290 * t289;
	t286 = sin(pkin(7));
	t287 = sin(pkin(6));
	t318 = t287 * t286;
	t276 = t288 * t313 - t318;
	t282 = t289 * pkin(9) + qJ(2);
	t293 = sin(qJ(3));
	t285 = sin(pkin(12));
	t297 = cos(qJ(3));
	t316 = t288 * t289;
	t296 = cos(qJ(4));
	t320 = t285 * t296;
	t292 = sin(qJ(4));
	t323 = pkin(11) * t292;
	t300 = (pkin(4) * t320 - pkin(10) * t316 + (pkin(3) + t323) * t285) * t297 - t288 * t286 * pkin(9) + t285 * pkin(2) + (pkin(3) * t316 + t285 * pkin(10)) * t293;
	t308 = pkin(4) * t296 + t323;
	t314 = t290 * t286;
	t274 = t287 * t289 + t288 * t314;
	t311 = t296 * t274;
	t312 = t292 * t274;
	t317 = t287 * t297;
	t329 = -t286 * pkin(10) * t317 + pkin(4) * t312 - pkin(11) * t311 + t287 * t282 + (pkin(3) * t318 - t308 * t276) * t293 - t300 * t290;
	t319 = t285 * t297;
	t303 = t276 * t293 + t290 * t319;
	t327 = t303 * t292 + t311;
	t324 = pkin(10) * t297;
	t275 = t287 * t316 + t314;
	t322 = t275 * t293;
	t321 = t285 * t293;
	t315 = t289 * t293;
	t310 = t297 * t288;
	t307 = -pkin(4) * t292 + pkin(11) * t296;
	t264 = t303 * t296 - t312;
	t268 = t276 * t297 - t290 * t321;
	t291 = sin(qJ(5));
	t295 = cos(qJ(5));
	t306 = t264 * t291 + t268 * t295;
	t305 = pkin(3) + t308;
	t304 = t285 * t317 + t322;
	t298 = cos(qJ(1));
	t294 = sin(qJ(1));
	t277 = t293 * t288 + t289 * t319;
	t273 = t288 * t318 - t313;
	t270 = t296 * t310 - t285 * (-t286 * t292 + t296 * t315);
	t269 = (-t285 * t315 + t310) * t292 - t286 * t320;
	t267 = t275 * t297 - t287 * t321;
	t266 = -t270 * t291 + t295 * t277;
	t265 = -t292 * t273 + t304 * t296;
	t263 = pkin(1) + (pkin(10) * t293 + t305 * t297 + pkin(2)) * t288 + ((pkin(9) - t307) * t286 + (-t305 * t293 + t324) * t289) * t285;
	t1 = [(-t264 * t294 + t298 * t270) * t295 + (t268 * t294 + t298 * t277) * t291, t266 * t298 + t306 * t294, t269 * t298 - t327 * t294, t263 * t298 + t329 * t294 + 0; (t264 * t295 - t268 * t291) * t298 + (t270 * t295 + t291 * t277) * t294, t266 * t294 - t306 * t298, t269 * t294 + t327 * t298, t263 * t294 - t329 * t298 + 0; t265 * t295 - t267 * t291, -t265 * t291 - t267 * t295, t296 * t273 + t304 * t292, pkin(8) + 0 + t308 * t322 + (t282 + (pkin(3) * t293 - t324) * t286) * t290 + t307 * t273 + t300 * t287; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:54:11
	% EndTime: 2020-11-04 21:54:12
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (186->73), mult. (445->119), div. (0->0), fcn. (580->14), ass. (0->54)
	t359 = cos(pkin(12));
	t360 = cos(pkin(7));
	t361 = cos(pkin(6));
	t378 = t361 * t360;
	t357 = sin(pkin(7));
	t358 = sin(pkin(6));
	t381 = t358 * t357;
	t346 = t359 * t378 - t381;
	t356 = sin(pkin(12));
	t365 = sin(qJ(3));
	t369 = cos(qJ(3));
	t377 = t361 * t369;
	t340 = t346 * t365 + t356 * t377;
	t379 = t361 * t357;
	t344 = t358 * t360 + t359 * t379;
	t364 = sin(qJ(4));
	t368 = cos(qJ(4));
	t391 = t340 * t364 + t368 * t344;
	t348 = t359 * t357 * pkin(9) - t356 * pkin(2);
	t380 = t359 * t360;
	t349 = -t356 * pkin(3) + pkin(10) * t380;
	t350 = pkin(3) * t380 + t356 * pkin(10);
	t352 = t360 * pkin(9) + qJ(2);
	t390 = (pkin(3) * t381 - t350 * t361) * t365 - (pkin(10) * t381 - t349 * t361) * t369 + t348 * t361 + t358 * t352;
	t366 = sin(qJ(1));
	t370 = cos(qJ(1));
	t382 = t356 * t370;
	t336 = (t346 * t366 + t360 * t382) * t369 - t365 * (t366 * t361 * t356 - t370 * t359);
	t363 = sin(qJ(5));
	t389 = t336 * t363;
	t383 = t356 * t369;
	t384 = t356 * t365;
	t337 = (-t346 * t369 + t361 * t384) * t370 + t366 * (t365 * t359 + t360 * t383);
	t388 = t337 * t363;
	t375 = t369 * t359;
	t339 = -t357 * t377 + (-t360 * t375 + t384) * t358;
	t387 = t339 * t363;
	t386 = t356 * t357;
	t385 = t356 * t360;
	t371 = (t358 * t380 + t379) * t365 + t358 * t383;
	t367 = cos(qJ(5));
	t362 = -qJ(6) - pkin(11);
	t353 = t367 * pkin(5) + pkin(4);
	t347 = -t360 * t384 + t375;
	t343 = t359 * t381 - t378;
	t342 = t347 * t364 - t368 * t386;
	t338 = (t359 * pkin(3) + pkin(10) * t385) * t369 + (-pkin(3) * t385 + t359 * pkin(10)) * t365 + pkin(9) * t386 + t359 * pkin(2) + pkin(1);
	t335 = -t368 * t343 - t371 * t364;
	t334 = -t364 * t343 + t371 * t368;
	t333 = (-t340 * t366 + t370 * t347) * t368 + (t344 * t366 + t357 * t382) * t364;
	t332 = -t342 * t366 - t391 * t370;
	t331 = (t340 * t368 - t364 * t344) * t370 + (t347 * t368 + t364 * t386) * t366;
	t330 = -t342 * t370 + t391 * t366;
	t1 = [t333 * t367 + t389, -t333 * t363 + t336 * t367, -t330, pkin(5) * t389 + t330 * t362 + t333 * t353 + t338 * t370 + t390 * t366 + 0; t331 * t367 + t388, -t331 * t363 + t337 * t367, -t332, pkin(5) * t388 + t331 * t353 + t332 * t362 + t338 * t366 - t390 * t370 + 0; t334 * t367 + t387, -t334 * t363 + t339 * t367, -t335, t334 * t353 + t335 * t362 + pkin(5) * t387 + (-pkin(10) * t379 - t349 * t358) * t369 + (pkin(3) * t379 + t350 * t358) * t365 - t348 * t358 + t352 * t361 + 0 + pkin(8); 0, 0, 0, 1;];
	Tc_mdh = t1;
end