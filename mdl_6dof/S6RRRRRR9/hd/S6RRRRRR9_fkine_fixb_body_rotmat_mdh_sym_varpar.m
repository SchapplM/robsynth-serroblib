% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRR9 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:49
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:57
	% EndTime: 2020-11-04 22:48:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:57
	% EndTime: 2020-11-04 22:48:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t166 = cos(qJ(1));
	t165 = sin(qJ(1));
	t1 = [t166, -t165, 0, 0; t165, t166, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:57
	% EndTime: 2020-11-04 22:48:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t167 = sin(pkin(6));
	t170 = sin(qJ(1));
	t178 = t170 * t167;
	t169 = sin(qJ(2));
	t177 = t170 * t169;
	t171 = cos(qJ(2));
	t176 = t170 * t171;
	t172 = cos(qJ(1));
	t175 = t172 * t167;
	t174 = t172 * t169;
	t173 = t172 * t171;
	t168 = cos(pkin(6));
	t1 = [-t168 * t177 + t173, -t168 * t176 - t174, t178, t172 * pkin(1) + pkin(9) * t178 + 0; t168 * t174 + t176, t168 * t173 - t177, -t175, t170 * pkin(1) - pkin(9) * t175 + 0; t167 * t169, t167 * t171, t168, t168 * pkin(9) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:57
	% EndTime: 2020-11-04 22:48:57
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->32), mult. (104->58), div. (0->0), fcn. (138->10), ass. (0->29)
	t185 = sin(pkin(7));
	t188 = cos(pkin(6));
	t207 = t185 * t188;
	t193 = cos(qJ(2));
	t206 = t185 * t193;
	t186 = sin(pkin(6));
	t187 = cos(pkin(7));
	t205 = t186 * t187;
	t204 = t188 * t193;
	t189 = sin(qJ(3));
	t190 = sin(qJ(2));
	t203 = t189 * t190;
	t202 = t189 * t193;
	t192 = cos(qJ(3));
	t201 = t190 * t192;
	t191 = sin(qJ(1));
	t200 = t191 * t190;
	t199 = t193 * t192;
	t194 = cos(qJ(1));
	t198 = t194 * t190;
	t197 = t194 * t193;
	t196 = -pkin(2) * t190 + pkin(10) * t206;
	t180 = -t185 * t186 + t187 * t204;
	t195 = t180 * t189 + t188 * t201;
	t184 = t187 * pkin(10) + pkin(9);
	t182 = t185 * t190 * pkin(10) + pkin(2) * t193 + pkin(1);
	t181 = -t187 * t203 + t199;
	t179 = t186 * t184 + t196 * t188;
	t1 = [t194 * t181 - t195 * t191, (-t180 * t191 - t187 * t198) * t192 - (-t188 * t200 + t197) * t189, (t191 * t204 + t198) * t185 + t191 * t205, t179 * t191 + t182 * t194 + 0; t191 * t181 + t195 * t194, (t180 * t192 - t188 * t203) * t194 - t191 * (t187 * t201 + t202), -(t188 * t197 - t200) * t185 - t194 * t205, -t179 * t194 + t182 * t191 + 0; t189 * t207 + (t187 * t202 + t201) * t186, t192 * t207 + (t187 * t199 - t203) * t186, -t186 * t206 + t188 * t187, t184 * t188 - t196 * t186 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:57
	% EndTime: 2020-11-04 22:48:57
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (91->44), mult. (209->80), div. (0->0), fcn. (264->12), ass. (0->38)
	t223 = sin(pkin(7));
	t225 = cos(pkin(7));
	t228 = sin(qJ(3));
	t232 = cos(qJ(3));
	t238 = pkin(3) * t228 - pkin(11) * t232;
	t213 = t223 * pkin(10) - t238 * t225;
	t218 = t232 * pkin(3) + pkin(11) * t228 + pkin(2);
	t229 = sin(qJ(2));
	t233 = cos(qJ(2));
	t250 = -t213 * t233 + t218 * t229;
	t224 = sin(pkin(6));
	t247 = t223 * t224;
	t246 = t223 * t228;
	t245 = t223 * t229;
	t226 = cos(pkin(6));
	t244 = t226 * t233;
	t243 = t228 * t229;
	t242 = t228 * t233;
	t241 = t229 * t232;
	t234 = cos(qJ(1));
	t240 = t229 * t234;
	t239 = t233 * t232;
	t236 = t225 * t242 + t241;
	t211 = -t224 * t246 + t236 * t226;
	t214 = t223 * t244 + t224 * t225;
	t227 = sin(qJ(4));
	t231 = cos(qJ(4));
	t237 = t211 * t227 + t231 * t214;
	t235 = t225 * pkin(10) + t238 * t223 + pkin(9);
	t230 = sin(qJ(1));
	t217 = t225 * t243 - t239;
	t216 = t225 * t244 - t247;
	t215 = -t226 * t225 + t233 * t247;
	t212 = t217 * t227 + t231 * t245;
	t210 = -t236 * t224 - t226 * t246;
	t209 = t213 * t229 + t218 * t233 + pkin(1);
	t208 = t224 * t235 - t250 * t226;
	t1 = [(-t211 * t230 - t234 * t217) * t231 + (t214 * t230 + t223 * t240) * t227, t234 * t212 + t237 * t230, (t216 * t230 + t225 * t240) * t232 + (-t230 * t226 * t229 + t234 * t233) * t228, t208 * t230 + t209 * t234 + 0; (t211 * t231 - t227 * t214) * t234 + t230 * (-t217 * t231 + t227 * t245), t230 * t212 - t237 * t234, (-t216 * t232 + t226 * t243) * t234 + t230 * (t225 * t241 + t242), -t208 * t234 + t209 * t230 + 0; -t210 * t231 - t227 * t215, t210 * t227 - t231 * t215, -t226 * t223 * t232 + (-t225 * t239 + t243) * t224, t250 * t224 + t235 * t226 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:57
	% EndTime: 2020-11-04 22:48:57
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (167->55), mult. (361->96), div. (0->0), fcn. (450->14), ass. (0->48)
	t271 = sin(pkin(7));
	t277 = sin(qJ(3));
	t282 = cos(qJ(3));
	t276 = sin(qJ(4));
	t281 = cos(qJ(4));
	t307 = t281 * pkin(4) + t276 * pkin(12) + pkin(3);
	t308 = -pkin(11) * t282 + t277 * t307;
	t310 = -t308 * t271 - pkin(9);
	t273 = cos(pkin(7));
	t289 = t276 * pkin(4) - t281 * pkin(12) + pkin(10);
	t256 = t271 * t289 - t273 * t308;
	t264 = pkin(11) * t277 + t282 * t307 + pkin(2);
	t278 = sin(qJ(2));
	t283 = cos(qJ(2));
	t309 = t256 * t283 - t264 * t278;
	t303 = t271 * t277;
	t302 = t271 * t282;
	t301 = t271 * t283;
	t272 = sin(pkin(6));
	t300 = t272 * t273;
	t299 = t277 * t278;
	t298 = t277 * t281;
	t297 = t277 * t283;
	t296 = t278 * t282;
	t295 = t283 * t282;
	t274 = cos(pkin(6));
	t287 = t271 * t298 + t273 * t276;
	t266 = -t271 * t276 + t273 * t298;
	t288 = t266 * t283 + t281 * t296;
	t254 = -t272 * t287 + t288 * t274;
	t286 = t273 * t295 - t299;
	t258 = -t272 * t302 + t286 * t274;
	t275 = sin(qJ(5));
	t280 = cos(qJ(5));
	t291 = t254 * t275 + t280 * t258;
	t285 = t273 * t297 + t296;
	t290 = (-t272 * t303 + t285 * t274) * t276 + t281 * (t274 * t301 + t300);
	t284 = cos(qJ(1));
	t279 = sin(qJ(1));
	t267 = t273 * t296 + t297;
	t261 = (t273 * t299 - t295) * t276 + t271 * t281 * t278;
	t260 = -t278 * t266 + t281 * t295;
	t259 = t286 * t272 + t274 * t302;
	t255 = -t260 * t275 + t280 * t267;
	t253 = -t288 * t272 - t274 * t287;
	t252 = t256 * t278 + t264 * t283 + pkin(1);
	t251 = -t310 * t272 + t309 * t274 + t289 * t300;
	t1 = [(-t254 * t279 + t284 * t260) * t280 + (t258 * t279 + t284 * t267) * t275, t284 * t255 + t291 * t279, -t284 * t261 - t290 * t279, t251 * t279 + t252 * t284 + 0; (t254 * t280 - t275 * t258) * t284 + t279 * (t260 * t280 + t275 * t267), t279 * t255 - t291 * t284, -t279 * t261 + t290 * t284, -t251 * t284 + t252 * t279 + 0; -t253 * t280 - t275 * t259, t253 * t275 - t280 * t259, (t285 * t272 + t274 * t303) * t276 + t281 * (t272 * t301 - t274 * t273), pkin(8) + 0 - t309 * t272 + (t289 * t273 - t310) * t274; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:57
	% EndTime: 2020-11-04 22:48:57
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (234->60), mult. (401->99), div. (0->0), fcn. (490->16), ass. (0->54)
	t330 = cos(qJ(5)) * pkin(5) + pkin(4);
	t338 = sin(qJ(4));
	t342 = cos(qJ(4));
	t346 = pkin(12) + pkin(13);
	t326 = t330 * t342 + t346 * t338 + pkin(3);
	t329 = sin(qJ(5)) * pkin(5) + pkin(11);
	t339 = sin(qJ(3));
	t343 = cos(qJ(3));
	t321 = t326 * t343 + t329 * t339 + pkin(2);
	t340 = sin(qJ(2));
	t366 = t321 * t340;
	t334 = sin(pkin(7));
	t365 = t334 * t339;
	t364 = t334 * t343;
	t344 = cos(qJ(2));
	t363 = t334 * t344;
	t362 = t339 * t342;
	t361 = t339 * t344;
	t360 = t340 * t339;
	t359 = t340 * t343;
	t358 = t344 * t343;
	t335 = sin(pkin(6));
	t337 = cos(pkin(6));
	t336 = cos(pkin(7));
	t351 = t334 * t362 + t336 * t338;
	t327 = -t334 * t338 + t336 * t362;
	t352 = t327 * t344 + t342 * t359;
	t315 = -t335 * t351 + t352 * t337;
	t322 = t340 * t327 - t342 * t358;
	t341 = sin(qJ(1));
	t345 = cos(qJ(1));
	t357 = t315 * t345 - t341 * t322;
	t356 = t315 * t341 + t345 * t322;
	t349 = t336 * t361 + t359;
	t355 = (-t335 * t365 + t349 * t337) * t338 + t342 * (t335 * t336 + t337 * t363);
	t354 = t326 * t339 - t329 * t343;
	t353 = t330 * t338 - t346 * t342 + pkin(10);
	t350 = t336 * t358 - t360;
	t348 = t354 * t336;
	t347 = t353 * t335;
	t333 = qJ(5) + qJ(6);
	t332 = cos(t333);
	t331 = sin(t333);
	t328 = t336 * t359 + t361;
	t323 = (t336 * t360 - t358) * t338 + t334 * t342 * t340;
	t320 = t350 * t335 + t337 * t364;
	t319 = -t335 * t364 + t350 * t337;
	t317 = t319 * t345 - t341 * t328;
	t316 = t319 * t341 + t345 * t328;
	t314 = t352 * t335 + t337 * t351;
	t313 = -t353 * t334 + t348;
	t312 = -t313 * t340 + t321 * t344 + pkin(1);
	t311 = (-t354 * t334 - pkin(9)) * t335 + (t313 * t344 + t366) * t337 - t336 * t347;
	t1 = [t316 * t331 - t356 * t332, t316 * t332 + t356 * t331, -t345 * t323 - t355 * t341, -t311 * t341 + t312 * t345 + 0; -t331 * t317 + t357 * t332, -t332 * t317 - t357 * t331, -t341 * t323 + t355 * t345, t311 * t345 + t312 * t341 + 0; t314 * t332 - t331 * t320, -t314 * t331 - t332 * t320, (t349 * t335 + t337 * t365) * t338 + t342 * (t335 * t363 - t337 * t336), pkin(8) + 0 + (t353 * t336 + pkin(9)) * t337 + (t344 * t348 + t366) * t335 + (t354 * t337 - t344 * t347) * t334; 0, 0, 0, 1;];
	Tc_mdh = t1;
end