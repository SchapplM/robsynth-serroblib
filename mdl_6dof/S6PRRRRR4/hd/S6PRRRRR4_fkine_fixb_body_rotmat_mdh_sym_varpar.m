% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRRR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:21
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRRRR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:33
	% EndTime: 2020-11-04 21:21:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:33
	% EndTime: 2020-11-04 21:21:33
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t148 = cos(pkin(13));
	t147 = sin(pkin(13));
	t1 = [t148, -t147, 0, 0; t147, t148, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:33
	% EndTime: 2020-11-04 21:21:33
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t149 = sin(pkin(13));
	t150 = sin(pkin(6));
	t158 = t149 * t150;
	t151 = cos(pkin(13));
	t157 = t151 * t150;
	t152 = cos(pkin(6));
	t153 = sin(qJ(2));
	t156 = t152 * t153;
	t154 = cos(qJ(2));
	t155 = t152 * t154;
	t1 = [-t149 * t156 + t151 * t154, -t149 * t155 - t151 * t153, t158, t151 * pkin(1) + pkin(8) * t158 + 0; t149 * t154 + t151 * t156, -t149 * t153 + t151 * t155, -t157, t149 * pkin(1) - pkin(8) * t157 + 0; t150 * t153, t150 * t154, t152, t152 * pkin(8) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:33
	% EndTime: 2020-11-04 21:21:33
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->35), mult. (112->66), div. (0->0), fcn. (146->10), ass. (0->30)
	t165 = sin(pkin(7));
	t187 = pkin(9) * t165;
	t164 = sin(pkin(13));
	t186 = t164 * pkin(2);
	t167 = cos(pkin(13));
	t185 = t167 * pkin(2);
	t169 = cos(pkin(6));
	t184 = t165 * t169;
	t173 = cos(qJ(2));
	t183 = t165 * t173;
	t168 = cos(pkin(7));
	t163 = t168 * pkin(9) + pkin(8);
	t166 = sin(pkin(6));
	t182 = t166 * t163;
	t181 = t166 * t168;
	t171 = sin(qJ(2));
	t180 = t168 * t171;
	t179 = t168 * t173;
	t178 = t169 * t171;
	t177 = t169 * t173;
	t176 = t164 * t187;
	t175 = t167 * t187;
	t174 = -t165 * t166 + t168 * t177;
	t172 = cos(qJ(3));
	t170 = sin(qJ(3));
	t162 = t164 * t178 - t167 * t173;
	t161 = t164 * t173 + t167 * t178;
	t160 = -t164 * t180 + t174 * t167;
	t159 = -t174 * t164 - t167 * t180;
	t1 = [t159 * t170 - t172 * t162, t159 * t172 + t170 * t162, (t164 * t177 + t167 * t171) * t165 + t164 * t181, (t169 * t176 + t185) * t173 + (-t169 * t186 + t175) * t171 + t164 * t182 + t167 * pkin(1) + 0; t160 * t170 + t161 * t172, t160 * t172 - t161 * t170, -(-t164 * t171 + t167 * t177) * t165 - t167 * t181, (-t169 * t175 + t186) * t173 + (t169 * t185 + t176) * t171 - t167 * t182 + t164 * pkin(1) + 0; t170 * t184 + (t170 * t179 + t171 * t172) * t166, t172 * t184 + (-t170 * t171 + t172 * t179) * t166, -t166 * t183 + t169 * t168, t163 * t169 + qJ(1) + 0 + (pkin(2) * t171 - pkin(9) * t183) * t166; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:33
	% EndTime: 2020-11-04 21:21:33
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (91->59), mult. (247->109), div. (0->0), fcn. (302->12), ass. (0->42)
	t199 = sin(pkin(6));
	t198 = sin(pkin(7));
	t201 = cos(pkin(7));
	t204 = sin(qJ(3));
	t207 = cos(qJ(3));
	t216 = pkin(3) * t204 - pkin(10) * t207;
	t213 = t201 * pkin(9) + t216 * t198 + pkin(8);
	t234 = t213 * t199;
	t205 = sin(qJ(2));
	t208 = cos(qJ(2));
	t232 = t198 * pkin(9);
	t212 = -t216 * t201 + t232;
	t215 = pkin(3) * t207 + pkin(10) * t204 + pkin(2);
	t233 = t215 * t205 - t212 * t208;
	t197 = sin(pkin(13));
	t231 = t197 * t201;
	t230 = t198 * t199;
	t229 = t199 * t201;
	t200 = cos(pkin(13));
	t202 = cos(pkin(6));
	t228 = t200 * t202;
	t227 = t201 * t204;
	t226 = t201 * t205;
	t225 = t201 * t208;
	t224 = t202 * t201;
	t223 = t202 * t204;
	t222 = t202 * t205;
	t221 = t202 * t207;
	t220 = t202 * t208;
	t219 = t200 * t224;
	t218 = t201 * t223;
	t217 = t204 * t230;
	t214 = t201 * t220 - t230;
	t211 = -(t197 * t218 - t200 * t207) * t208 - (t197 * t221 + t200 * t227) * t205 + t197 * t217;
	t210 = -(t197 * t207 + t200 * t218) * t208 - (-t197 * t227 + t200 * t221) * t205 + t200 * t217;
	t206 = cos(qJ(4));
	t203 = sin(qJ(4));
	t195 = t208 * t230 - t224;
	t190 = t198 * t223 + (t204 * t225 + t207 * t205) * t199;
	t189 = t200 * t229 + (-t197 * t205 + t200 * t220) * t198;
	t188 = t198 * t205 * t200 + (t198 * t220 + t229) * t197;
	t1 = [t203 * t188 + t211 * t206, t206 * t188 - t211 * t203, (t214 * t197 + t200 * t226) * t207 - t204 * (t197 * t222 - t200 * t208), 0 + (t212 * t205 + t215 * t208 + pkin(1)) * t200 + (-t233 * t202 + t234) * t197; -t203 * t189 - t210 * t206, -t206 * t189 + t210 * t203, (t197 * t226 - t214 * t200) * t207 + (t197 * t208 + t200 * t222) * t204, ((t197 * pkin(3) - pkin(10) * t219) * t207 + (pkin(3) * t219 + t197 * pkin(10)) * t204 - t228 * t232 + t197 * pkin(2)) * t208 + ((pkin(3) * t228 + pkin(10) * t231) * t207 + (-pkin(3) * t231 + pkin(10) * t228) * t204 + pkin(2) * t228 + t197 * t232) * t205 + t197 * pkin(1) + 0 - t200 * t234; t190 * t206 - t203 * t195, -t190 * t203 - t206 * t195, -t198 * t221 + (t204 * t205 - t207 * t225) * t199, t233 * t199 + t213 * t202 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:33
	% EndTime: 2020-11-04 21:21:33
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (117->55), mult. (262->89), div. (0->0), fcn. (342->14), ass. (0->45)
	t279 = pkin(4) * sin(qJ(4));
	t254 = sin(pkin(7));
	t278 = pkin(9) * t254;
	t253 = sin(pkin(13));
	t277 = t253 * pkin(2);
	t256 = cos(pkin(13));
	t276 = t256 * pkin(2);
	t258 = cos(pkin(6));
	t275 = t254 * t258;
	t263 = cos(qJ(2));
	t274 = t254 * t263;
	t257 = cos(pkin(7));
	t248 = t257 * pkin(9) + pkin(8);
	t255 = sin(pkin(6));
	t273 = t255 * t248;
	t272 = t255 * t257;
	t261 = sin(qJ(2));
	t271 = t257 * t261;
	t270 = t257 * t263;
	t269 = t258 * t261;
	t268 = t258 * t263;
	t267 = t253 * t278;
	t266 = t256 * t278;
	t265 = -t254 * t255 + t257 * t268;
	t264 = -pkin(11) - pkin(10);
	t262 = cos(qJ(3));
	t260 = sin(qJ(3));
	t252 = qJ(4) + qJ(5);
	t251 = cos(t252);
	t250 = sin(t252);
	t249 = cos(qJ(4)) * pkin(4) + pkin(3);
	t247 = t253 * t269 - t256 * t263;
	t246 = t253 * t263 + t256 * t269;
	t245 = -t255 * t274 + t258 * t257;
	t244 = (t253 * t268 + t256 * t261) * t254 + t253 * t272;
	t243 = -(-t253 * t261 + t256 * t268) * t254 - t256 * t272;
	t242 = t260 * t275 + (t260 * t270 + t261 * t262) * t255;
	t241 = -t262 * t275 + (t260 * t261 - t262 * t270) * t255;
	t240 = -t253 * t271 + t265 * t256;
	t239 = -t265 * t253 - t256 * t271;
	t238 = t240 * t262 - t246 * t260;
	t237 = t240 * t260 + t246 * t262;
	t236 = t239 * t262 + t260 * t247;
	t235 = t239 * t260 - t262 * t247;
	t1 = [t235 * t251 + t244 * t250, -t235 * t250 + t244 * t251, -t236, t235 * t249 + t236 * t264 + t244 * t279 + (t258 * t267 + t276) * t263 + (-t258 * t277 + t266) * t261 + t253 * t273 + t256 * pkin(1) + 0; t237 * t251 + t243 * t250, -t237 * t250 + t243 * t251, -t238, t237 * t249 + t238 * t264 + t243 * t279 + (-t258 * t266 + t277) * t263 + (t258 * t276 + t267) * t261 - t256 * t273 + t253 * pkin(1) + 0; t242 * t251 + t245 * t250, -t242 * t250 + t245 * t251, t241, t245 * t279 - t241 * t264 + t242 * t249 + t248 * t258 + qJ(1) + 0 + (pkin(2) * t261 - pkin(9) * t274) * t255; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:33
	% EndTime: 2020-11-04 21:21:33
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (211->65), mult. (467->107), div. (0->0), fcn. (617->16), ass. (0->53)
	t333 = pkin(4) * sin(qJ(4));
	t305 = sin(pkin(7));
	t332 = pkin(9) * t305;
	t304 = sin(pkin(13));
	t331 = t304 * pkin(2);
	t307 = cos(pkin(13));
	t330 = t307 * pkin(2);
	t309 = cos(pkin(6));
	t329 = t305 * t309;
	t316 = cos(qJ(2));
	t328 = t305 * t316;
	t308 = cos(pkin(7));
	t299 = t308 * pkin(9) + pkin(8);
	t306 = sin(pkin(6));
	t327 = t306 * t299;
	t326 = t306 * t308;
	t313 = sin(qJ(2));
	t325 = t308 * t313;
	t324 = t308 * t316;
	t323 = t309 * t313;
	t322 = t309 * t316;
	t321 = t304 * t332;
	t320 = t307 * t332;
	t319 = t304 * t316 + t307 * t323;
	t318 = -t305 * t306 + t308 * t322;
	t317 = -pkin(11) - pkin(10);
	t315 = cos(qJ(3));
	t314 = cos(qJ(6));
	t312 = sin(qJ(3));
	t310 = sin(qJ(6));
	t303 = qJ(4) + qJ(5);
	t302 = cos(t303);
	t301 = sin(t303);
	t300 = cos(qJ(4)) * pkin(4) + pkin(3);
	t298 = t304 * t323 - t307 * t316;
	t297 = -t306 * t328 + t309 * t308;
	t295 = (t304 * t322 + t307 * t313) * t305 + t304 * t326;
	t294 = -(-t304 * t313 + t307 * t322) * t305 - t307 * t326;
	t293 = t312 * t329 + (t312 * t324 + t313 * t315) * t306;
	t292 = -t315 * t329 + (t312 * t313 - t315 * t324) * t306;
	t291 = -t304 * t325 + t318 * t307;
	t290 = -t318 * t304 - t307 * t325;
	t289 = t291 * t315 - t319 * t312;
	t288 = t291 * t312 + t319 * t315;
	t287 = t290 * t315 + t312 * t298;
	t286 = t290 * t312 - t315 * t298;
	t285 = t293 * t302 + t297 * t301;
	t284 = t293 * t301 - t297 * t302;
	t283 = t288 * t302 + t294 * t301;
	t282 = t288 * t301 - t294 * t302;
	t281 = t286 * t302 + t295 * t301;
	t280 = t286 * t301 - t295 * t302;
	t1 = [t281 * t314 - t287 * t310, -t281 * t310 - t287 * t314, t280, t281 * pkin(5) + t280 * pkin(12) + t286 * t300 + t287 * t317 + t295 * t333 + (t309 * t321 + t330) * t316 + (-t309 * t331 + t320) * t313 + t304 * t327 + t307 * pkin(1) + 0; t283 * t314 - t289 * t310, -t283 * t310 - t289 * t314, t282, t283 * pkin(5) + t282 * pkin(12) + t288 * t300 + t289 * t317 + t294 * t333 + (-t309 * t320 + t331) * t316 + (t309 * t330 + t321) * t313 - t307 * t327 + t304 * pkin(1) + 0; t285 * t314 + t292 * t310, -t285 * t310 + t292 * t314, t284, t297 * t333 + t285 * pkin(5) + t284 * pkin(12) - t292 * t317 + t293 * t300 + t299 * t309 + qJ(1) + 0 + (pkin(2) * t313 - pkin(9) * t328) * t306; 0, 0, 0, 1;];
	Tc_mdh = t1;
end