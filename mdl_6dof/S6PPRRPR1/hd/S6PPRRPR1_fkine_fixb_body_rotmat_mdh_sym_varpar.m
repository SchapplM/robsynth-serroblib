% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PPRRPR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:53
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PPRRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:53:44
	% EndTime: 2020-11-04 20:53:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:53:44
	% EndTime: 2020-11-04 20:53:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t151 = cos(pkin(11));
	t150 = sin(pkin(11));
	t1 = [t151, -t150, 0, 0; t150, t151, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:53:44
	% EndTime: 2020-11-04 20:53:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t153 = sin(pkin(11));
	t154 = sin(pkin(6));
	t161 = t153 * t154;
	t157 = cos(pkin(6));
	t160 = t153 * t157;
	t156 = cos(pkin(11));
	t159 = t156 * t154;
	t158 = t156 * t157;
	t155 = cos(pkin(12));
	t152 = sin(pkin(12));
	t1 = [-t152 * t160 + t156 * t155, -t156 * t152 - t155 * t160, t161, t156 * pkin(1) + qJ(2) * t161 + 0; t152 * t158 + t153 * t155, -t153 * t152 + t155 * t158, -t159, t153 * pkin(1) - qJ(2) * t159 + 0; t154 * t152, t154 * t155, t157, t157 * qJ(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:53:44
	% EndTime: 2020-11-04 20:53:44
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (44->30), mult. (101->54), div. (0->0), fcn. (135->10), ass. (0->29)
	t173 = sin(pkin(7));
	t190 = pkin(8) * t173;
	t171 = sin(pkin(12));
	t176 = cos(pkin(11));
	t189 = t171 * t176;
	t172 = sin(pkin(11));
	t188 = t172 * t171;
	t178 = cos(pkin(6));
	t187 = t173 * t178;
	t174 = sin(pkin(6));
	t186 = t174 * t173;
	t175 = cos(pkin(12));
	t177 = cos(pkin(7));
	t185 = t175 * t177;
	t184 = t176 * t178;
	t183 = t177 * t174;
	t182 = t178 * t177;
	t168 = -t171 * pkin(2) + t175 * t190;
	t169 = pkin(8) * t177 + qJ(2);
	t181 = t168 * t178 + t169 * t174;
	t180 = cos(qJ(3));
	t179 = sin(qJ(3));
	t167 = pkin(2) * t175 + t171 * t190 + pkin(1);
	t166 = t175 * t176 - t178 * t188;
	t165 = t171 * t184 + t172 * t175;
	t164 = t175 * t182 - t186;
	t163 = -t164 * t172 - t177 * t189;
	t162 = t164 * t176 - t177 * t188;
	t1 = [t163 * t179 + t166 * t180, t163 * t180 - t166 * t179, (t175 * t187 + t183) * t172 + t173 * t189, t167 * t176 + t172 * t181 + 0; t162 * t179 + t165 * t180, t162 * t180 - t165 * t179, (-t175 * t184 + t188) * t173 - t176 * t183, t167 * t172 - t176 * t181 + 0; t179 * t187 + (t171 * t180 + t179 * t185) * t174, t180 * t187 + (-t171 * t179 + t180 * t185) * t174, -t175 * t186 + t182, -t168 * t174 + t169 * t178 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:53:44
	% EndTime: 2020-11-04 20:53:45
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (91->53), mult. (217->95), div. (0->0), fcn. (272->12), ass. (0->43)
	t211 = sin(pkin(12));
	t213 = sin(pkin(7));
	t236 = t211 * t213;
	t217 = cos(pkin(7));
	t235 = t211 * t217;
	t218 = cos(pkin(6));
	t234 = t211 * t218;
	t214 = sin(pkin(6));
	t233 = t214 * t213;
	t215 = cos(pkin(12));
	t232 = t215 * t217;
	t231 = t217 * t214;
	t230 = t218 * t213;
	t229 = t218 * t217;
	t197 = t215 * t229 - t233;
	t212 = sin(pkin(11));
	t216 = cos(pkin(11));
	t193 = t197 * t212 + t216 * t235;
	t200 = -t212 * t234 + t216 * t215;
	t220 = sin(qJ(3));
	t222 = cos(qJ(3));
	t228 = t193 * t220 - t200 * t222;
	t194 = -t197 * t216 + t212 * t235;
	t199 = t212 * t215 + t216 * t234;
	t227 = t194 * t220 - t199 * t222;
	t202 = t215 * t213 * pkin(8) - t211 * pkin(2);
	t208 = t217 * pkin(8) + qJ(2);
	t226 = t202 * t218 + t214 * t208;
	t206 = pkin(3) * t232 + t211 * pkin(9);
	t225 = pkin(3) * t233 - t206 * t218;
	t205 = -t211 * pkin(3) + pkin(9) * t232;
	t224 = pkin(9) * t233 - t205 * t218;
	t223 = t214 * t211 * t222 + (t215 * t231 + t230) * t220;
	t221 = cos(qJ(4));
	t219 = sin(qJ(4));
	t204 = pkin(3) * t235 - t215 * pkin(9);
	t203 = t215 * pkin(3) + pkin(9) * t235;
	t201 = t215 * pkin(2) + pkin(8) * t236 + pkin(1);
	t196 = t215 * t233 - t229;
	t195 = t215 * t230 + t231;
	t192 = t195 * t216 - t212 * t236;
	t191 = t195 * t212 + t216 * t236;
	t1 = [t191 * t219 - t228 * t221, t191 * t221 + t228 * t219, t193 * t222 + t200 * t220, (t216 * t203 - t224 * t212) * t222 + (-t216 * t204 + t225 * t212) * t220 + t226 * t212 + t201 * t216 + 0; -t192 * t219 - t227 * t221, -t192 * t221 + t227 * t219, t194 * t222 + t199 * t220, (t212 * t203 + t224 * t216) * t222 + (-t212 * t204 - t225 * t216) * t220 - t226 * t216 + t201 * t212 + 0; -t219 * t196 + t223 * t221, -t221 * t196 - t223 * t219, -t222 * t230 + (t211 * t220 - t222 * t232) * t214, (-pkin(9) * t230 - t205 * t214) * t222 + (pkin(3) * t230 + t206 * t214) * t220 - t202 * t214 + t208 * t218 + 0 + qJ(1); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:53:45
	% EndTime: 2020-11-04 20:53:45
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (167->68), mult. (411->113), div. (0->0), fcn. (536->14), ass. (0->54)
	t267 = sin(pkin(12));
	t269 = sin(pkin(7));
	t293 = t267 * t269;
	t274 = cos(pkin(7));
	t292 = t267 * t274;
	t275 = cos(pkin(6));
	t291 = t267 * t275;
	t270 = sin(pkin(6));
	t290 = t270 * t269;
	t272 = cos(pkin(12));
	t289 = t272 * t274;
	t288 = t274 * t270;
	t287 = t275 * t269;
	t286 = t275 * t274;
	t252 = t272 * t286 - t290;
	t268 = sin(pkin(11));
	t273 = cos(pkin(11));
	t248 = t252 * t268 + t273 * t292;
	t255 = -t268 * t291 + t273 * t272;
	t277 = sin(qJ(3));
	t279 = cos(qJ(3));
	t285 = t248 * t277 - t255 * t279;
	t249 = -t252 * t273 + t268 * t292;
	t254 = t268 * t272 + t273 * t291;
	t284 = t249 * t277 - t254 * t279;
	t257 = t272 * t269 * pkin(8) - t267 * pkin(2);
	t263 = t274 * pkin(8) + qJ(2);
	t283 = t257 * t275 + t270 * t263;
	t261 = pkin(3) * t289 + t267 * pkin(9);
	t282 = pkin(3) * t290 - t261 * t275;
	t260 = -t267 * pkin(3) + pkin(9) * t289;
	t281 = pkin(9) * t290 - t260 * t275;
	t280 = t270 * t267 * t279 + (t272 * t288 + t287) * t277;
	t278 = cos(qJ(4));
	t276 = sin(qJ(4));
	t271 = cos(pkin(13));
	t266 = sin(pkin(13));
	t259 = pkin(3) * t292 - t272 * pkin(9);
	t258 = t272 * pkin(3) + pkin(9) * t292;
	t256 = t272 * pkin(2) + pkin(8) * t293 + pkin(1);
	t251 = t272 * t290 - t286;
	t250 = t272 * t287 + t288;
	t247 = t250 * t273 - t268 * t293;
	t246 = t250 * t268 + t273 * t293;
	t245 = -t279 * t287 + (t267 * t277 - t279 * t289) * t270;
	t244 = -t278 * t251 - t280 * t276;
	t243 = -t276 * t251 + t280 * t278;
	t242 = t249 * t279 + t254 * t277;
	t241 = t248 * t279 + t255 * t277;
	t240 = -t247 * t278 + t284 * t276;
	t239 = t246 * t276 - t285 * t278;
	t238 = -t247 * t276 - t284 * t278;
	t237 = t246 * t278 + t285 * t276;
	t1 = [t239 * t271 + t241 * t266, -t239 * t266 + t241 * t271, -t237, t239 * pkin(4) - t237 * qJ(5) + (t273 * t258 - t281 * t268) * t279 + (-t273 * t259 + t282 * t268) * t277 + t283 * t268 + t256 * t273 + 0; t238 * t271 + t242 * t266, -t238 * t266 + t242 * t271, -t240, t238 * pkin(4) - t240 * qJ(5) + (t268 * t258 + t281 * t273) * t279 + (-t268 * t259 - t282 * t273) * t277 - t283 * t273 + t256 * t268 + 0; t243 * t271 + t245 * t266, -t243 * t266 + t245 * t271, -t244, t243 * pkin(4) - t244 * qJ(5) + (-pkin(9) * t287 - t260 * t270) * t279 + (pkin(3) * t287 + t261 * t270) * t277 - t257 * t270 + t263 * t275 + 0 + qJ(1); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:53:45
	% EndTime: 2020-11-04 20:53:45
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (210->75), mult. (443->124), div. (0->0), fcn. (556->16), ass. (0->60)
	t319 = cos(pkin(13)) * pkin(5) + pkin(4);
	t333 = qJ(5) + pkin(10);
	t334 = sin(qJ(4));
	t336 = cos(qJ(4));
	t357 = -t319 * t334 + t333 * t336;
	t325 = sin(pkin(12));
	t327 = sin(pkin(7));
	t355 = t325 * t327;
	t328 = sin(pkin(6));
	t354 = t325 * t328;
	t331 = cos(pkin(7));
	t353 = t325 * t331;
	t332 = cos(pkin(6));
	t352 = t325 * t332;
	t351 = t328 * t327;
	t329 = cos(pkin(12));
	t350 = t329 * t331;
	t349 = t331 * t328;
	t348 = t332 * t327;
	t347 = t332 * t331;
	t337 = cos(qJ(3));
	t345 = t337 * t354;
	t311 = t329 * t347 - t351;
	t326 = sin(pkin(11));
	t330 = cos(pkin(11));
	t301 = t311 * t326 + t330 * t353;
	t308 = t326 * t352 - t329 * t330;
	t335 = sin(qJ(3));
	t344 = t301 * t335 + t308 * t337;
	t302 = t311 * t330 - t326 * t353;
	t313 = t326 * t329 + t330 * t352;
	t343 = t302 * t335 + t313 * t337;
	t315 = pkin(8) * t327 * t329 - t325 * pkin(2);
	t317 = pkin(8) * t331 + qJ(2);
	t342 = t315 * t332 + t328 * t317;
	t341 = t319 * t336 + t333 * t334;
	t318 = sin(pkin(13)) * pkin(5) + pkin(9);
	t305 = pkin(3) * t350 + t318 * t325;
	t340 = pkin(3) * t351 - t305 * t332;
	t304 = -t325 * pkin(3) + t318 * t350;
	t339 = -t304 * t332 + t318 * t351;
	t312 = t329 * t349 + t348;
	t338 = t312 * t335 + t345;
	t324 = pkin(13) + qJ(6);
	t321 = cos(t324);
	t320 = sin(t324);
	t314 = pkin(2) * t329 + pkin(8) * t355 + pkin(1);
	t310 = t329 * t351 - t347;
	t309 = t329 * t348 + t349;
	t307 = pkin(3) * t353 - t318 * t329;
	t306 = pkin(3) * t329 + t318 * t353;
	t303 = t312 * t337 - t335 * t354;
	t300 = t309 * t330 - t326 * t355;
	t299 = t309 * t326 + t330 * t355;
	t298 = t302 * t337 - t313 * t335;
	t297 = t301 * t337 - t308 * t335;
	t296 = -t310 * t334 + t336 * t338;
	t295 = -t300 * t334 + t336 * t343;
	t294 = -t299 * t334 + t336 * t344;
	t1 = [-t294 * t321 + t297 * t320, t294 * t320 + t297 * t321, -t299 * t336 - t334 * t344, (-t301 * t341 - t307 * t330 + t326 * t340) * t335 + (t330 * t306 - t308 * t341 - t326 * t339) * t337 + t342 * t326 + t314 * t330 + 0 - t357 * t299; t295 * t321 - t298 * t320, -t295 * t320 - t298 * t321, t300 * t336 + t334 * t343, (t302 * t341 - t307 * t326 - t330 * t340) * t335 + (t306 * t326 + t313 * t341 + t330 * t339) * t337 - t342 * t330 + t314 * t326 + 0 + t357 * t300; t296 * t321 - t303 * t320, -t296 * t320 - t303 * t321, t336 * t310 + t334 * t338, (pkin(3) * t348 + t305 * t328 + t312 * t341) * t335 + (t310 * t333 + t319 * t345) * t336 + (-t310 * t319 + t333 * t345) * t334 + (-t304 * t328 - t318 * t348) * t337 - t315 * t328 + t317 * t332 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end