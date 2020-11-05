% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRP11 (for one body)
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
% Datum: 2020-11-04 22:45
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRRP11_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:45:43
	% EndTime: 2020-11-04 22:45:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:45:43
	% EndTime: 2020-11-04 22:45:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t159 = cos(qJ(1));
	t158 = sin(qJ(1));
	t1 = [t159, -t158, 0, 0; t158, t159, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:45:43
	% EndTime: 2020-11-04 22:45:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t160 = sin(pkin(6));
	t163 = sin(qJ(1));
	t171 = t163 * t160;
	t162 = sin(qJ(2));
	t170 = t163 * t162;
	t164 = cos(qJ(2));
	t169 = t163 * t164;
	t165 = cos(qJ(1));
	t168 = t165 * t160;
	t167 = t165 * t162;
	t166 = t165 * t164;
	t161 = cos(pkin(6));
	t1 = [-t161 * t170 + t166, -t161 * t169 - t167, t171, t165 * pkin(1) + pkin(9) * t171 + 0; t161 * t167 + t169, t161 * t166 - t170, -t168, t163 * pkin(1) - pkin(9) * t168 + 0; t160 * t162, t160 * t164, t161, t161 * pkin(9) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:45:43
	% EndTime: 2020-11-04 22:45:43
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (44->32), mult. (104->58), div. (0->0), fcn. (138->10), ass. (0->29)
	t178 = sin(pkin(7));
	t181 = cos(pkin(6));
	t200 = t178 * t181;
	t186 = cos(qJ(2));
	t199 = t178 * t186;
	t179 = sin(pkin(6));
	t180 = cos(pkin(7));
	t198 = t179 * t180;
	t197 = t181 * t186;
	t182 = sin(qJ(3));
	t183 = sin(qJ(2));
	t196 = t182 * t183;
	t195 = t182 * t186;
	t185 = cos(qJ(3));
	t194 = t183 * t185;
	t184 = sin(qJ(1));
	t193 = t184 * t183;
	t192 = t185 * t186;
	t187 = cos(qJ(1));
	t191 = t187 * t183;
	t190 = t187 * t186;
	t189 = -pkin(2) * t183 + pkin(10) * t199;
	t173 = -t178 * t179 + t180 * t197;
	t188 = t173 * t182 + t181 * t194;
	t177 = t180 * pkin(10) + pkin(9);
	t175 = t178 * t183 * pkin(10) + pkin(2) * t186 + pkin(1);
	t174 = -t180 * t196 + t192;
	t172 = t179 * t177 + t189 * t181;
	t1 = [t187 * t174 - t188 * t184, (-t173 * t184 - t180 * t191) * t185 - (-t181 * t193 + t190) * t182, (t184 * t197 + t191) * t178 + t184 * t198, t172 * t184 + t175 * t187 + 0; t184 * t174 + t188 * t187, (t173 * t185 - t181 * t196) * t187 - t184 * (t180 * t194 + t195), -(t181 * t190 - t193) * t178 - t187 * t198, -t172 * t187 + t175 * t184 + 0; t182 * t200 + (t180 * t195 + t194) * t179, t185 * t200 + (t180 * t192 - t196) * t179, -t179 * t199 + t181 * t180, t177 * t181 - t189 * t179 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:45:43
	% EndTime: 2020-11-04 22:45:44
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (91->44), mult. (209->80), div. (0->0), fcn. (264->12), ass. (0->38)
	t216 = sin(pkin(7));
	t218 = cos(pkin(7));
	t221 = sin(qJ(3));
	t225 = cos(qJ(3));
	t231 = pkin(3) * t221 - pkin(11) * t225;
	t206 = t216 * pkin(10) - t231 * t218;
	t211 = t225 * pkin(3) + pkin(11) * t221 + pkin(2);
	t222 = sin(qJ(2));
	t226 = cos(qJ(2));
	t243 = -t206 * t226 + t211 * t222;
	t217 = sin(pkin(6));
	t240 = t216 * t217;
	t239 = t216 * t221;
	t238 = t216 * t222;
	t219 = cos(pkin(6));
	t237 = t219 * t226;
	t236 = t221 * t222;
	t235 = t221 * t226;
	t234 = t222 * t225;
	t227 = cos(qJ(1));
	t233 = t222 * t227;
	t232 = t225 * t226;
	t229 = t218 * t235 + t234;
	t203 = -t217 * t239 + t229 * t219;
	t207 = t216 * t237 + t217 * t218;
	t220 = sin(qJ(4));
	t224 = cos(qJ(4));
	t230 = t203 * t220 + t224 * t207;
	t228 = t218 * pkin(10) + t231 * t216 + pkin(9);
	t223 = sin(qJ(1));
	t210 = t218 * t236 - t232;
	t209 = t218 * t237 - t240;
	t208 = -t219 * t218 + t226 * t240;
	t205 = t210 * t220 + t224 * t238;
	t204 = t229 * t217 + t219 * t239;
	t202 = t206 * t222 + t211 * t226 + pkin(1);
	t201 = t217 * t228 - t243 * t219;
	t1 = [(-t203 * t223 - t227 * t210) * t224 + (t207 * t223 + t216 * t233) * t220, t227 * t205 + t230 * t223, (t209 * t223 + t218 * t233) * t225 + (-t223 * t219 * t222 + t227 * t226) * t221, t201 * t223 + t202 * t227 + 0; (t203 * t224 - t220 * t207) * t227 + (-t210 * t224 + t220 * t238) * t223, t205 * t223 - t230 * t227, (-t209 * t225 + t219 * t236) * t227 + t223 * (t218 * t234 + t235), -t201 * t227 + t202 * t223 + 0; t204 * t224 - t220 * t208, -t204 * t220 - t224 * t208, -t219 * t216 * t225 + (-t218 * t232 + t236) * t217, t243 * t217 + t228 * t219 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:45:44
	% EndTime: 2020-11-04 22:45:44
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (167->55), mult. (361->96), div. (0->0), fcn. (450->14), ass. (0->48)
	t264 = sin(pkin(7));
	t270 = sin(qJ(3));
	t275 = cos(qJ(3));
	t269 = sin(qJ(4));
	t274 = cos(qJ(4));
	t300 = t274 * pkin(4) + t269 * pkin(12) + pkin(3);
	t301 = -t275 * pkin(11) + t270 * t300;
	t303 = -t301 * t264 - pkin(9);
	t266 = cos(pkin(7));
	t282 = t269 * pkin(4) - t274 * pkin(12) + pkin(10);
	t249 = t264 * t282 - t266 * t301;
	t257 = pkin(11) * t270 + t275 * t300 + pkin(2);
	t271 = sin(qJ(2));
	t276 = cos(qJ(2));
	t302 = t249 * t276 - t257 * t271;
	t296 = t264 * t270;
	t295 = t264 * t275;
	t294 = t264 * t276;
	t265 = sin(pkin(6));
	t293 = t265 * t266;
	t292 = t270 * t271;
	t291 = t270 * t274;
	t290 = t270 * t276;
	t289 = t271 * t275;
	t288 = t275 * t276;
	t267 = cos(pkin(6));
	t280 = t264 * t291 + t266 * t269;
	t259 = -t264 * t269 + t266 * t291;
	t281 = t259 * t276 + t274 * t289;
	t247 = -t265 * t280 + t281 * t267;
	t279 = t266 * t288 - t292;
	t251 = -t265 * t295 + t279 * t267;
	t268 = sin(qJ(5));
	t273 = cos(qJ(5));
	t284 = t247 * t268 + t273 * t251;
	t278 = t266 * t290 + t289;
	t283 = (-t265 * t296 + t278 * t267) * t269 + t274 * (t267 * t294 + t293);
	t277 = cos(qJ(1));
	t272 = sin(qJ(1));
	t260 = t266 * t289 + t290;
	t254 = (t266 * t292 - t288) * t269 + t264 * t274 * t271;
	t253 = -t271 * t259 + t274 * t288;
	t252 = t279 * t265 + t267 * t295;
	t248 = -t253 * t268 + t273 * t260;
	t246 = -t281 * t265 - t267 * t280;
	t245 = t249 * t271 + t257 * t276 + pkin(1);
	t244 = -t303 * t265 + t302 * t267 + t282 * t293;
	t1 = [(-t247 * t272 + t253 * t277) * t273 + t268 * (t251 * t272 + t277 * t260), t277 * t248 + t284 * t272, -t277 * t254 - t283 * t272, t244 * t272 + t245 * t277 + 0; (t247 * t273 - t268 * t251) * t277 + t272 * (t253 * t273 + t268 * t260), t272 * t248 - t284 * t277, -t254 * t272 + t283 * t277, -t244 * t277 + t245 * t272 + 0; -t246 * t273 - t268 * t252, t246 * t268 - t273 * t252, (t278 * t265 + t267 * t296) * t269 + t274 * (t265 * t294 - t267 * t266), pkin(8) + 0 - t302 * t265 + (t282 * t266 - t303) * t267; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:45:44
	% EndTime: 2020-11-04 22:45:44
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (222->61), mult. (401->103), div. (0->0), fcn. (490->14), ass. (0->51)
	t334 = cos(qJ(5));
	t323 = t334 * pkin(5) + pkin(4);
	t328 = pkin(12) + qJ(6);
	t330 = sin(qJ(4));
	t335 = cos(qJ(4));
	t318 = t323 * t335 + t328 * t330 + pkin(3);
	t329 = sin(qJ(5));
	t322 = t329 * pkin(5) + pkin(11);
	t331 = sin(qJ(3));
	t336 = cos(qJ(3));
	t313 = t318 * t336 + t322 * t331 + pkin(2);
	t332 = sin(qJ(2));
	t357 = t313 * t332;
	t324 = sin(pkin(7));
	t356 = t324 * t331;
	t355 = t324 * t336;
	t337 = cos(qJ(2));
	t354 = t324 * t337;
	t353 = t331 * t335;
	t352 = t331 * t337;
	t351 = t332 * t331;
	t350 = t332 * t336;
	t349 = t336 * t337;
	t325 = sin(pkin(6));
	t327 = cos(pkin(6));
	t326 = cos(pkin(7));
	t343 = t324 * t353 + t326 * t330;
	t320 = -t324 * t330 + t326 * t353;
	t344 = t320 * t337 + t335 * t350;
	t308 = -t325 * t343 + t344 * t327;
	t342 = t326 * t349 - t351;
	t311 = -t325 * t355 + t342 * t327;
	t348 = t308 * t329 + t334 * t311;
	t341 = t326 * t352 + t350;
	t347 = (-t325 * t356 + t341 * t327) * t330 + t335 * (t325 * t326 + t327 * t354);
	t346 = t318 * t331 - t322 * t336;
	t345 = t323 * t330 - t328 * t335 + pkin(10);
	t340 = t346 * t326;
	t339 = t325 * t345;
	t338 = cos(qJ(1));
	t333 = sin(qJ(1));
	t321 = t326 * t350 + t352;
	t315 = (t326 * t351 - t349) * t330 + t324 * t335 * t332;
	t314 = -t332 * t320 + t335 * t349;
	t312 = t342 * t325 + t327 * t355;
	t309 = -t314 * t329 + t334 * t321;
	t307 = -t344 * t325 - t327 * t343;
	t306 = -t324 * t345 + t340;
	t305 = -t306 * t332 + t313 * t337 + pkin(1);
	t304 = t325 * (-t346 * t324 - pkin(9)) + (t306 * t337 + t357) * t327 - t326 * t339;
	t1 = [(-t308 * t333 + t314 * t338) * t334 + t329 * (t311 * t333 + t338 * t321), t338 * t309 + t348 * t333, -t338 * t315 - t347 * t333, -t304 * t333 + t305 * t338 + 0; (t308 * t334 - t329 * t311) * t338 + t333 * (t314 * t334 + t329 * t321), t333 * t309 - t348 * t338, -t315 * t333 + t347 * t338, t304 * t338 + t305 * t333 + 0; -t307 * t334 - t329 * t312, t307 * t329 - t334 * t312, (t341 * t325 + t327 * t356) * t330 + t335 * (t325 * t354 - t327 * t326), 0 + pkin(8) + (t345 * t326 + pkin(9)) * t327 + (t337 * t340 + t357) * t325 + (t346 * t327 - t337 * t339) * t324; 0, 0, 0, 1;];
	Tc_mdh = t1;
end