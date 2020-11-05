% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRP12 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:46
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRRP12_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP12_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:46:09
	% EndTime: 2020-11-04 22:46:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:46:09
	% EndTime: 2020-11-04 22:46:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t163 = cos(qJ(1));
	t162 = sin(qJ(1));
	t1 = [t163, -t162, 0, 0; t162, t163, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:46:09
	% EndTime: 2020-11-04 22:46:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t164 = sin(pkin(6));
	t167 = sin(qJ(1));
	t175 = t167 * t164;
	t166 = sin(qJ(2));
	t174 = t167 * t166;
	t168 = cos(qJ(2));
	t173 = t167 * t168;
	t169 = cos(qJ(1));
	t172 = t169 * t164;
	t171 = t169 * t166;
	t170 = t169 * t168;
	t165 = cos(pkin(6));
	t1 = [-t165 * t174 + t170, -t165 * t173 - t171, t175, t169 * pkin(1) + pkin(9) * t175 + 0; t165 * t171 + t173, t165 * t170 - t174, -t172, t167 * pkin(1) - pkin(9) * t172 + 0; t164 * t166, t164 * t168, t165, t165 * pkin(9) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:46:09
	% EndTime: 2020-11-04 22:46:09
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (44->32), mult. (104->58), div. (0->0), fcn. (138->10), ass. (0->29)
	t182 = sin(pkin(7));
	t185 = cos(pkin(6));
	t204 = t182 * t185;
	t190 = cos(qJ(2));
	t203 = t182 * t190;
	t183 = sin(pkin(6));
	t184 = cos(pkin(7));
	t202 = t183 * t184;
	t201 = t185 * t190;
	t186 = sin(qJ(3));
	t187 = sin(qJ(2));
	t200 = t186 * t187;
	t199 = t186 * t190;
	t189 = cos(qJ(3));
	t198 = t187 * t189;
	t188 = sin(qJ(1));
	t197 = t188 * t187;
	t196 = t189 * t190;
	t191 = cos(qJ(1));
	t195 = t191 * t187;
	t194 = t191 * t190;
	t193 = -pkin(2) * t187 + pkin(10) * t203;
	t177 = -t182 * t183 + t184 * t201;
	t192 = t177 * t186 + t185 * t198;
	t181 = t184 * pkin(10) + pkin(9);
	t179 = t182 * t187 * pkin(10) + pkin(2) * t190 + pkin(1);
	t178 = -t184 * t200 + t196;
	t176 = t183 * t181 + t193 * t185;
	t1 = [t191 * t178 - t192 * t188, (-t177 * t188 - t184 * t195) * t189 - (-t185 * t197 + t194) * t186, (t188 * t201 + t195) * t182 + t188 * t202, t176 * t188 + t179 * t191 + 0; t188 * t178 + t192 * t191, (t177 * t189 - t185 * t200) * t191 - t188 * (t184 * t198 + t199), -(t185 * t194 - t197) * t182 - t191 * t202, -t176 * t191 + t179 * t188 + 0; t186 * t204 + (t184 * t199 + t198) * t183, t189 * t204 + (t184 * t196 - t200) * t183, -t183 * t203 + t185 * t184, t181 * t185 - t193 * t183 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:46:09
	% EndTime: 2020-11-04 22:46:09
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (91->44), mult. (209->80), div. (0->0), fcn. (264->12), ass. (0->38)
	t220 = sin(pkin(7));
	t222 = cos(pkin(7));
	t225 = sin(qJ(3));
	t229 = cos(qJ(3));
	t235 = pkin(3) * t225 - pkin(11) * t229;
	t210 = t220 * pkin(10) - t235 * t222;
	t215 = t229 * pkin(3) + pkin(11) * t225 + pkin(2);
	t226 = sin(qJ(2));
	t230 = cos(qJ(2));
	t247 = -t210 * t230 + t215 * t226;
	t221 = sin(pkin(6));
	t244 = t220 * t221;
	t243 = t220 * t225;
	t242 = t220 * t226;
	t223 = cos(pkin(6));
	t241 = t223 * t230;
	t240 = t225 * t226;
	t239 = t225 * t230;
	t238 = t226 * t229;
	t231 = cos(qJ(1));
	t237 = t226 * t231;
	t236 = t229 * t230;
	t233 = t222 * t239 + t238;
	t208 = -t221 * t243 + t233 * t223;
	t211 = t220 * t241 + t221 * t222;
	t224 = sin(qJ(4));
	t228 = cos(qJ(4));
	t234 = t208 * t224 + t228 * t211;
	t232 = t222 * pkin(10) + t235 * t220 + pkin(9);
	t227 = sin(qJ(1));
	t214 = t222 * t240 - t236;
	t213 = t222 * t241 - t244;
	t212 = -t223 * t222 + t230 * t244;
	t209 = t214 * t224 + t228 * t242;
	t207 = -t233 * t221 - t223 * t243;
	t206 = t210 * t226 + t215 * t230 + pkin(1);
	t205 = t221 * t232 - t247 * t223;
	t1 = [(-t208 * t227 - t231 * t214) * t228 + (t211 * t227 + t220 * t237) * t224, t209 * t231 + t234 * t227, (t213 * t227 + t222 * t237) * t229 + (-t227 * t223 * t226 + t231 * t230) * t225, t205 * t227 + t206 * t231 + 0; (t208 * t228 - t224 * t211) * t231 + t227 * (-t214 * t228 + t224 * t242), t209 * t227 - t234 * t231, (-t213 * t229 + t223 * t240) * t231 + t227 * (t222 * t238 + t239), -t205 * t231 + t206 * t227 + 0; -t207 * t228 - t224 * t212, t207 * t224 - t228 * t212, -t223 * t220 * t229 + (-t222 * t236 + t240) * t221, t247 * t221 + t232 * t223 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:46:09
	% EndTime: 2020-11-04 22:46:10
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (167->55), mult. (361->96), div. (0->0), fcn. (450->14), ass. (0->48)
	t268 = sin(pkin(7));
	t274 = sin(qJ(3));
	t279 = cos(qJ(3));
	t273 = sin(qJ(4));
	t278 = cos(qJ(4));
	t304 = t278 * pkin(4) + t273 * pkin(12) + pkin(3);
	t305 = -pkin(11) * t279 + t274 * t304;
	t307 = -t305 * t268 - pkin(9);
	t270 = cos(pkin(7));
	t286 = t273 * pkin(4) - t278 * pkin(12) + pkin(10);
	t253 = t268 * t286 - t270 * t305;
	t261 = pkin(11) * t274 + t279 * t304 + pkin(2);
	t275 = sin(qJ(2));
	t280 = cos(qJ(2));
	t306 = t253 * t280 - t261 * t275;
	t300 = t268 * t274;
	t299 = t268 * t279;
	t298 = t268 * t280;
	t269 = sin(pkin(6));
	t297 = t269 * t270;
	t296 = t274 * t275;
	t295 = t274 * t278;
	t294 = t274 * t280;
	t293 = t275 * t279;
	t292 = t279 * t280;
	t271 = cos(pkin(6));
	t284 = t268 * t295 + t273 * t270;
	t263 = -t273 * t268 + t270 * t295;
	t285 = t263 * t280 + t278 * t293;
	t251 = -t269 * t284 + t285 * t271;
	t283 = t270 * t292 - t296;
	t255 = -t269 * t299 + t283 * t271;
	t272 = sin(qJ(5));
	t277 = cos(qJ(5));
	t288 = t251 * t272 + t277 * t255;
	t282 = t270 * t294 + t293;
	t287 = (-t269 * t300 + t282 * t271) * t273 + t278 * (t271 * t298 + t297);
	t281 = cos(qJ(1));
	t276 = sin(qJ(1));
	t264 = t270 * t293 + t294;
	t258 = (t270 * t296 - t292) * t273 + t268 * t278 * t275;
	t257 = -t275 * t263 + t278 * t292;
	t256 = t283 * t269 + t271 * t299;
	t252 = -t257 * t272 + t277 * t264;
	t250 = -t285 * t269 - t271 * t284;
	t249 = t253 * t275 + t261 * t280 + pkin(1);
	t248 = -t307 * t269 + t306 * t271 + t286 * t297;
	t1 = [(-t251 * t276 + t281 * t257) * t277 + t272 * (t255 * t276 + t281 * t264), t252 * t281 + t288 * t276, -t258 * t281 - t287 * t276, t248 * t276 + t249 * t281 + 0; (t251 * t277 - t272 * t255) * t281 + t276 * (t257 * t277 + t272 * t264), t276 * t252 - t288 * t281, -t258 * t276 + t287 * t281, -t248 * t281 + t249 * t276 + 0; -t250 * t277 - t272 * t256, t250 * t272 - t277 * t256, (t282 * t269 + t271 * t300) * t273 + t278 * (t269 * t298 - t271 * t270), pkin(8) + 0 - t306 * t269 + (t286 * t270 - t307) * t271; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:46:10
	% EndTime: 2020-11-04 22:46:10
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (235->59), mult. (429->100), div. (0->0), fcn. (518->14), ass. (0->50)
	t330 = sin(pkin(7));
	t334 = sin(qJ(5));
	t339 = cos(qJ(5));
	t327 = -t334 * pkin(5) + qJ(6) * t339 - pkin(11);
	t336 = sin(qJ(3));
	t341 = cos(qJ(3));
	t326 = pkin(5) * t339 + qJ(6) * t334 + pkin(4);
	t335 = sin(qJ(4));
	t340 = cos(qJ(4));
	t366 = t335 * pkin(12) + t326 * t340 + pkin(3);
	t349 = t327 * t341 + t336 * t366;
	t368 = -t349 * t330 - pkin(9);
	t332 = cos(pkin(7));
	t348 = t340 * pkin(12) - t326 * t335 - pkin(10);
	t310 = t330 * t348 + t349 * t332;
	t314 = -t327 * t336 + t341 * t366 + pkin(2);
	t337 = sin(qJ(2));
	t342 = cos(qJ(2));
	t367 = t310 * t342 + t314 * t337;
	t362 = t330 * t336;
	t361 = t330 * t341;
	t360 = t330 * t342;
	t331 = sin(pkin(6));
	t359 = t331 * t332;
	t358 = t336 * t337;
	t357 = t336 * t340;
	t356 = t336 * t342;
	t355 = t337 * t341;
	t354 = t341 * t342;
	t333 = cos(pkin(6));
	t346 = t330 * t357 + t335 * t332;
	t324 = -t335 * t330 + t332 * t357;
	t347 = t324 * t342 + t340 * t355;
	t312 = -t331 * t346 + t347 * t333;
	t345 = t332 * t354 - t358;
	t316 = -t331 * t361 + t345 * t333;
	t351 = t312 * t334 + t339 * t316;
	t344 = t332 * t356 + t355;
	t350 = (-t331 * t362 + t344 * t333) * t335 + t340 * (t333 * t360 + t359);
	t343 = cos(qJ(1));
	t338 = sin(qJ(1));
	t325 = t332 * t355 + t356;
	t319 = (t332 * t358 - t354) * t335 + t330 * t340 * t337;
	t318 = -t337 * t324 + t340 * t354;
	t317 = t345 * t331 + t333 * t361;
	t313 = -t318 * t334 + t339 * t325;
	t311 = t347 * t331 + t333 * t346;
	t309 = -t310 * t337 + t314 * t342 + pkin(1);
	t308 = t331 * t368 + t367 * t333 + t348 * t359;
	t1 = [(-t312 * t338 + t343 * t318) * t339 + t334 * (t316 * t338 + t343 * t325), -t319 * t343 - t350 * t338, -t313 * t343 - t351 * t338, -t308 * t338 + t309 * t343 + 0; (t312 * t339 - t334 * t316) * t343 + t338 * (t318 * t339 + t334 * t325), -t319 * t338 + t350 * t343, -t338 * t313 + t351 * t343, t308 * t343 + t309 * t338 + 0; t311 * t339 - t334 * t317, (t344 * t331 + t333 * t362) * t335 + t340 * (t331 * t360 - t333 * t332), t311 * t334 + t339 * t317, pkin(8) + 0 + t367 * t331 + (-t348 * t332 - t368) * t333; 0, 0, 0, 1;];
	Tc_mdh = t1;
end