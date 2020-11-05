% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PPRRRR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:55
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PPRRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:55:43
	% EndTime: 2020-11-04 20:55:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:55:43
	% EndTime: 2020-11-04 20:55:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t137 = cos(pkin(12));
	t136 = sin(pkin(12));
	t1 = [t137, -t136, 0, 0; t136, t137, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:55:43
	% EndTime: 2020-11-04 20:55:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t139 = sin(pkin(12));
	t140 = sin(pkin(6));
	t147 = t139 * t140;
	t143 = cos(pkin(6));
	t146 = t139 * t143;
	t142 = cos(pkin(12));
	t145 = t142 * t140;
	t144 = t142 * t143;
	t141 = cos(pkin(13));
	t138 = sin(pkin(13));
	t1 = [-t138 * t146 + t142 * t141, -t142 * t138 - t141 * t146, t147, t142 * pkin(1) + qJ(2) * t147 + 0; t138 * t144 + t139 * t141, -t139 * t138 + t141 * t144, -t145, t139 * pkin(1) - qJ(2) * t145 + 0; t140 * t138, t140 * t141, t143, t143 * qJ(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:55:43
	% EndTime: 2020-11-04 20:55:44
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (44->30), mult. (101->54), div. (0->0), fcn. (135->10), ass. (0->29)
	t159 = sin(pkin(7));
	t176 = pkin(8) * t159;
	t157 = sin(pkin(13));
	t162 = cos(pkin(12));
	t175 = t157 * t162;
	t158 = sin(pkin(12));
	t174 = t158 * t157;
	t164 = cos(pkin(6));
	t173 = t159 * t164;
	t160 = sin(pkin(6));
	t172 = t160 * t159;
	t161 = cos(pkin(13));
	t163 = cos(pkin(7));
	t171 = t161 * t163;
	t170 = t162 * t164;
	t169 = t163 * t160;
	t168 = t164 * t163;
	t154 = -t157 * pkin(2) + t161 * t176;
	t155 = t163 * pkin(8) + qJ(2);
	t167 = t154 * t164 + t160 * t155;
	t166 = cos(qJ(3));
	t165 = sin(qJ(3));
	t153 = t161 * pkin(2) + t157 * t176 + pkin(1);
	t152 = t162 * t161 - t164 * t174;
	t151 = t157 * t170 + t158 * t161;
	t150 = t161 * t168 - t172;
	t149 = -t150 * t158 - t163 * t175;
	t148 = t150 * t162 - t163 * t174;
	t1 = [t149 * t165 + t152 * t166, t149 * t166 - t152 * t165, (t161 * t173 + t169) * t158 + t159 * t175, t153 * t162 + t167 * t158 + 0; t148 * t165 + t151 * t166, t148 * t166 - t151 * t165, (-t161 * t170 + t174) * t159 - t162 * t169, t153 * t158 - t167 * t162 + 0; t165 * t173 + (t157 * t166 + t165 * t171) * t160, t166 * t173 + (-t157 * t165 + t166 * t171) * t160, -t161 * t172 + t168, -t154 * t160 + t155 * t164 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:55:44
	% EndTime: 2020-11-04 20:55:44
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (91->53), mult. (217->95), div. (0->0), fcn. (272->12), ass. (0->43)
	t197 = sin(pkin(13));
	t199 = sin(pkin(7));
	t222 = t197 * t199;
	t203 = cos(pkin(7));
	t221 = t197 * t203;
	t204 = cos(pkin(6));
	t220 = t197 * t204;
	t200 = sin(pkin(6));
	t219 = t200 * t199;
	t201 = cos(pkin(13));
	t218 = t201 * t203;
	t217 = t203 * t200;
	t216 = t204 * t199;
	t215 = t204 * t203;
	t184 = t201 * t215 - t219;
	t198 = sin(pkin(12));
	t202 = cos(pkin(12));
	t179 = t184 * t198 + t202 * t221;
	t186 = -t198 * t220 + t202 * t201;
	t206 = sin(qJ(3));
	t208 = cos(qJ(3));
	t214 = t179 * t206 - t186 * t208;
	t180 = -t184 * t202 + t198 * t221;
	t185 = t198 * t201 + t202 * t220;
	t213 = t180 * t206 - t185 * t208;
	t188 = t201 * t199 * pkin(8) - t197 * pkin(2);
	t194 = t203 * pkin(8) + qJ(2);
	t212 = t188 * t204 + t200 * t194;
	t192 = pkin(3) * t218 + t197 * pkin(9);
	t211 = pkin(3) * t219 - t192 * t204;
	t191 = -t197 * pkin(3) + pkin(9) * t218;
	t210 = pkin(9) * t219 - t191 * t204;
	t209 = t200 * t197 * t208 + (t201 * t217 + t216) * t206;
	t207 = cos(qJ(4));
	t205 = sin(qJ(4));
	t190 = pkin(3) * t221 - t201 * pkin(9);
	t189 = t201 * pkin(3) + pkin(9) * t221;
	t187 = t201 * pkin(2) + pkin(8) * t222 + pkin(1);
	t182 = t201 * t219 - t215;
	t181 = t201 * t216 + t217;
	t178 = t181 * t202 - t198 * t222;
	t177 = t181 * t198 + t202 * t222;
	t1 = [t177 * t205 - t214 * t207, t177 * t207 + t214 * t205, t179 * t208 + t186 * t206, (t202 * t189 - t210 * t198) * t208 + (-t202 * t190 + t211 * t198) * t206 + t212 * t198 + t187 * t202 + 0; -t205 * t178 - t213 * t207, -t207 * t178 + t213 * t205, t180 * t208 + t185 * t206, (t198 * t189 + t210 * t202) * t208 + (-t198 * t190 - t211 * t202) * t206 - t212 * t202 + t187 * t198 + 0; -t205 * t182 + t209 * t207, -t207 * t182 - t209 * t205, -t208 * t216 + (t197 * t206 - t208 * t218) * t200, (-pkin(9) * t216 - t191 * t200) * t208 + (pkin(3) * t216 + t192 * t200) * t206 - t188 * t200 + t194 * t204 + 0 + qJ(1); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:55:44
	% EndTime: 2020-11-04 20:55:44
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (117->56), mult. (242->89), div. (0->0), fcn. (322->14), ass. (0->50)
	t275 = pkin(4) * sin(qJ(4));
	t249 = sin(pkin(7));
	t274 = pkin(8) * t249;
	t247 = sin(pkin(13));
	t251 = cos(pkin(13));
	t252 = cos(pkin(12));
	t248 = sin(pkin(12));
	t254 = cos(pkin(6));
	t269 = t248 * t254;
	t231 = t247 * t269 - t252 * t251;
	t257 = cos(qJ(3));
	t273 = t231 * t257;
	t266 = t252 * t254;
	t236 = t247 * t266 + t248 * t251;
	t256 = sin(qJ(3));
	t272 = t236 * t256;
	t271 = t247 * t256;
	t270 = t248 * t247;
	t250 = sin(pkin(6));
	t268 = t250 * t249;
	t267 = t252 * t247;
	t253 = cos(pkin(7));
	t265 = t253 * t250;
	t264 = t254 * t249;
	t263 = t254 * t253;
	t262 = t256 * t231;
	t235 = t251 * t263 - t268;
	t227 = -t235 * t252 + t253 * t270;
	t230 = t236 * t257;
	t261 = t227 * t256 - t230;
	t234 = t251 * t265 + t264;
	t260 = t250 * t247 * t257 + t234 * t256;
	t226 = t235 * t248 + t253 * t267;
	t259 = t226 * t256 + t273;
	t258 = pkin(10) + pkin(9);
	t246 = qJ(4) + qJ(5);
	t245 = cos(t246);
	t244 = sin(t246);
	t243 = cos(qJ(4)) * pkin(4) + pkin(3);
	t242 = t253 * pkin(8) + qJ(2);
	t238 = t247 * pkin(2) - t251 * t274;
	t237 = t251 * pkin(2) + t247 * t274 + pkin(1);
	t233 = t251 * t268 - t263;
	t232 = t251 * t264 + t265;
	t229 = t238 * t254 - t242 * t250;
	t228 = (t251 * t266 - t270) * t253 - t252 * t268;
	t225 = t232 * t252 - t249 * t270;
	t224 = t232 * t248 + t249 * t267;
	t223 = (t251 * t269 + t267) * t253 - t248 * t268;
	t1 = [t224 * t244 - t259 * t245, t224 * t245 + t259 * t244, t226 * t257 - t262, -(t256 * t223 + t273) * t243 + (t223 * t257 - t262) * t258 + t224 * t275 - t229 * t248 + t237 * t252 + 0; -t225 * t244 - t261 * t245, -t245 * t225 + t261 * t244, t227 * t257 + t272, (t256 * t228 + t230) * t243 - (t228 * t257 - t272) * t258 - t225 * t275 + t252 * t229 + t237 * t248 + 0; -t244 * t233 + t260 * t245, -t245 * t233 - t260 * t244, -t257 * t264 + (-t251 * t253 * t257 + t271) * t250, t260 * t243 - t258 * (t234 * t257 - t250 * t271) - t233 * t275 + t242 * t254 + t238 * t250 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:55:44
	% EndTime: 2020-11-04 20:55:44
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (211->71), mult. (436->107), div. (0->0), fcn. (586->16), ass. (0->61)
	t339 = pkin(4) * sin(qJ(4));
	t311 = sin(pkin(7));
	t338 = pkin(8) * t311;
	t309 = sin(pkin(13));
	t313 = cos(pkin(13));
	t314 = cos(pkin(12));
	t310 = sin(pkin(12));
	t316 = cos(pkin(6));
	t333 = t310 * t316;
	t293 = t309 * t333 - t314 * t313;
	t321 = cos(qJ(3));
	t337 = t293 * t321;
	t330 = t314 * t316;
	t298 = t309 * t330 + t310 * t313;
	t319 = sin(qJ(3));
	t336 = t298 * t319;
	t335 = t309 * t319;
	t334 = t310 * t309;
	t312 = sin(pkin(6));
	t332 = t312 * t311;
	t331 = t314 * t309;
	t315 = cos(pkin(7));
	t329 = t315 * t312;
	t328 = t316 * t311;
	t327 = t316 * t315;
	t326 = t319 * t293;
	t297 = t313 * t327 - t332;
	t289 = -t297 * t314 + t315 * t334;
	t292 = t298 * t321;
	t325 = t289 * t319 - t292;
	t296 = t313 * t329 + t328;
	t324 = t312 * t309 * t321 + t296 * t319;
	t288 = t297 * t310 + t315 * t331;
	t323 = t288 * t319 + t337;
	t322 = pkin(10) + pkin(9);
	t320 = cos(qJ(6));
	t317 = sin(qJ(6));
	t308 = qJ(4) + qJ(5);
	t307 = cos(t308);
	t306 = sin(t308);
	t305 = cos(qJ(4)) * pkin(4) + pkin(3);
	t304 = t315 * pkin(8) + qJ(2);
	t300 = t309 * pkin(2) - t313 * t338;
	t299 = t313 * pkin(2) + t309 * t338 + pkin(1);
	t295 = t313 * t332 - t327;
	t294 = t313 * t328 + t329;
	t291 = t300 * t316 - t304 * t312;
	t290 = (t313 * t330 - t334) * t315 - t314 * t332;
	t287 = t294 * t314 - t311 * t334;
	t286 = t294 * t310 + t311 * t331;
	t285 = (t313 * t333 + t331) * t315 - t310 * t332;
	t284 = -t321 * t328 + (-t313 * t315 * t321 + t335) * t312;
	t283 = t289 * t321 + t336;
	t282 = t288 * t321 - t326;
	t281 = -t307 * t295 - t324 * t306;
	t280 = -t306 * t295 + t324 * t307;
	t279 = -t307 * t287 + t325 * t306;
	t278 = t286 * t306 - t323 * t307;
	t277 = -t287 * t306 - t325 * t307;
	t276 = t286 * t307 + t323 * t306;
	t1 = [t278 * t320 + t282 * t317, -t278 * t317 + t282 * t320, -t276, t278 * pkin(5) - t276 * pkin(11) - (t319 * t285 + t337) * t305 + (t285 * t321 - t326) * t322 + t286 * t339 - t291 * t310 + t299 * t314 + 0; t277 * t320 + t283 * t317, -t277 * t317 + t283 * t320, -t279, t277 * pkin(5) - t279 * pkin(11) + (t319 * t290 + t292) * t305 - (t290 * t321 - t336) * t322 - t287 * t339 + t314 * t291 + t299 * t310 + 0; t280 * t320 + t284 * t317, -t280 * t317 + t284 * t320, -t281, t280 * pkin(5) - t281 * pkin(11) + t324 * t305 - t322 * (t296 * t321 - t312 * t335) - t295 * t339 + t304 * t316 + t300 * t312 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end