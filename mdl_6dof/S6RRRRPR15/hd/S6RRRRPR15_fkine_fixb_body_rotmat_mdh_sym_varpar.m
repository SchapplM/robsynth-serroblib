% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPR15 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:42
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRPR15_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR15_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:42:17
	% EndTime: 2020-11-04 22:42:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:42:17
	% EndTime: 2020-11-04 22:42:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t158 = cos(qJ(1));
	t157 = sin(qJ(1));
	t1 = [t158, -t157, 0, 0; t157, t158, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:42:17
	% EndTime: 2020-11-04 22:42:17
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t159 = sin(pkin(6));
	t162 = sin(qJ(1));
	t170 = t162 * t159;
	t161 = sin(qJ(2));
	t169 = t162 * t161;
	t163 = cos(qJ(2));
	t168 = t162 * t163;
	t164 = cos(qJ(1));
	t167 = t164 * t159;
	t166 = t164 * t161;
	t165 = t164 * t163;
	t160 = cos(pkin(6));
	t1 = [-t160 * t169 + t165, -t160 * t168 - t166, t170, t164 * pkin(1) + pkin(9) * t170 + 0; t160 * t166 + t168, t160 * t165 - t169, -t167, t162 * pkin(1) - pkin(9) * t167 + 0; t159 * t161, t159 * t163, t160, t160 * pkin(9) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:42:17
	% EndTime: 2020-11-04 22:42:17
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->32), mult. (104->58), div. (0->0), fcn. (138->10), ass. (0->29)
	t177 = sin(pkin(7));
	t180 = cos(pkin(6));
	t199 = t177 * t180;
	t185 = cos(qJ(2));
	t198 = t177 * t185;
	t178 = sin(pkin(6));
	t179 = cos(pkin(7));
	t197 = t178 * t179;
	t196 = t180 * t185;
	t181 = sin(qJ(3));
	t182 = sin(qJ(2));
	t195 = t181 * t182;
	t194 = t181 * t185;
	t184 = cos(qJ(3));
	t193 = t182 * t184;
	t183 = sin(qJ(1));
	t192 = t183 * t182;
	t191 = t184 * t185;
	t186 = cos(qJ(1));
	t190 = t186 * t182;
	t189 = t186 * t185;
	t188 = -pkin(2) * t182 + pkin(10) * t198;
	t172 = -t177 * t178 + t179 * t196;
	t187 = t172 * t181 + t180 * t193;
	t176 = t179 * pkin(10) + pkin(9);
	t174 = t177 * t182 * pkin(10) + pkin(2) * t185 + pkin(1);
	t173 = -t179 * t195 + t191;
	t171 = t178 * t176 + t188 * t180;
	t1 = [t186 * t173 - t187 * t183, (-t172 * t183 - t179 * t190) * t184 - (-t180 * t192 + t189) * t181, (t183 * t196 + t190) * t177 + t183 * t197, t171 * t183 + t174 * t186 + 0; t183 * t173 + t187 * t186, (t172 * t184 - t180 * t195) * t186 - t183 * (t179 * t193 + t194), -(t180 * t189 - t192) * t177 - t186 * t197, -t171 * t186 + t174 * t183 + 0; t181 * t199 + (t179 * t194 + t193) * t178, t184 * t199 + (t179 * t191 - t195) * t178, -t178 * t198 + t180 * t179, t176 * t180 - t188 * t178 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:42:17
	% EndTime: 2020-11-04 22:42:17
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (91->44), mult. (209->80), div. (0->0), fcn. (264->12), ass. (0->38)
	t215 = sin(pkin(7));
	t217 = cos(pkin(7));
	t220 = sin(qJ(3));
	t224 = cos(qJ(3));
	t230 = pkin(3) * t220 - pkin(11) * t224;
	t205 = t215 * pkin(10) - t230 * t217;
	t210 = t224 * pkin(3) + pkin(11) * t220 + pkin(2);
	t221 = sin(qJ(2));
	t225 = cos(qJ(2));
	t242 = -t205 * t225 + t210 * t221;
	t216 = sin(pkin(6));
	t239 = t215 * t216;
	t238 = t215 * t220;
	t237 = t215 * t221;
	t218 = cos(pkin(6));
	t236 = t218 * t225;
	t235 = t220 * t221;
	t234 = t220 * t225;
	t233 = t221 * t224;
	t226 = cos(qJ(1));
	t232 = t221 * t226;
	t231 = t224 * t225;
	t228 = t217 * t234 + t233;
	t203 = -t216 * t238 + t228 * t218;
	t206 = t215 * t236 + t216 * t217;
	t219 = sin(qJ(4));
	t223 = cos(qJ(4));
	t229 = t203 * t219 + t223 * t206;
	t227 = t217 * pkin(10) + t230 * t215 + pkin(9);
	t222 = sin(qJ(1));
	t209 = t217 * t235 - t231;
	t208 = t217 * t236 - t239;
	t207 = -t218 * t217 + t225 * t239;
	t204 = t209 * t219 + t223 * t237;
	t202 = -t228 * t216 - t218 * t238;
	t201 = t205 * t221 + t210 * t225 + pkin(1);
	t200 = t216 * t227 - t242 * t218;
	t1 = [(-t203 * t222 - t226 * t209) * t223 + t219 * (t206 * t222 + t215 * t232), t204 * t226 + t229 * t222, (t208 * t222 + t217 * t232) * t224 + (-t222 * t218 * t221 + t226 * t225) * t220, t200 * t222 + t201 * t226 + 0; (t203 * t223 - t219 * t206) * t226 + t222 * (-t209 * t223 + t219 * t237), t204 * t222 - t229 * t226, (-t208 * t224 + t218 * t235) * t226 + t222 * (t217 * t233 + t234), -t200 * t226 + t201 * t222 + 0; -t202 * t223 - t219 * t207, t202 * t219 - t223 * t207, -t218 * t215 * t224 + (-t217 * t231 + t235) * t216, t242 * t216 + t227 * t218 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:42:17
	% EndTime: 2020-11-04 22:42:18
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (133->49), mult. (264->85), div. (0->0), fcn. (319->12), ass. (0->41)
	t259 = sin(pkin(7));
	t264 = sin(qJ(3));
	t268 = cos(qJ(3));
	t263 = sin(qJ(4));
	t267 = cos(qJ(4));
	t290 = pkin(4) * t267 + qJ(5) * t263 + pkin(3);
	t291 = -t268 * pkin(11) + t264 * t290;
	t293 = -t291 * t259 - pkin(9);
	t261 = cos(pkin(7));
	t272 = t263 * pkin(4) - qJ(5) * t267 + pkin(10);
	t245 = -t259 * t272 + t261 * t291;
	t249 = pkin(11) * t264 + t268 * t290 + pkin(2);
	t265 = sin(qJ(2));
	t269 = cos(qJ(2));
	t292 = t245 * t269 + t249 * t265;
	t260 = sin(pkin(6));
	t286 = t259 * t260;
	t285 = t259 * t264;
	t284 = t259 * t265;
	t283 = t260 * t261;
	t262 = cos(pkin(6));
	t282 = t262 * t269;
	t281 = t264 * t265;
	t280 = t264 * t269;
	t279 = t265 * t268;
	t270 = cos(qJ(1));
	t278 = t265 * t270;
	t277 = t268 * t269;
	t271 = t261 * t280 + t279;
	t247 = -t260 * t285 + t271 * t262;
	t250 = t259 * t282 + t283;
	t273 = t247 * t263 + t267 * t250;
	t266 = sin(qJ(1));
	t253 = t261 * t281 - t277;
	t252 = t261 * t282 - t286;
	t251 = -t262 * t261 + t269 * t286;
	t248 = t253 * t263 + t267 * t284;
	t246 = -t271 * t260 - t262 * t285;
	t244 = -t245 * t265 + t249 * t269 + pkin(1);
	t243 = t260 * t293 + t292 * t262 - t272 * t283;
	t1 = [(t252 * t266 + t261 * t278) * t268 + (-t266 * t262 * t265 + t270 * t269) * t264, (t247 * t266 + t270 * t253) * t267 - t263 * (t250 * t266 + t259 * t278), -t248 * t270 - t273 * t266, -t243 * t266 + t244 * t270 + 0; (-t252 * t268 + t262 * t281) * t270 + t266 * (t261 * t279 + t280), (-t247 * t267 + t263 * t250) * t270 - t266 * (-t253 * t267 + t263 * t284), -t248 * t266 + t273 * t270, t243 * t270 + t244 * t266 + 0; -t262 * t259 * t268 + (-t261 * t277 + t281) * t260, t246 * t267 + t263 * t251, -t246 * t263 + t267 * t251, pkin(8) + 0 + t292 * t260 + (t272 * t261 - t293) * t262; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:42:18
	% EndTime: 2020-11-04 22:42:18
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (201->58), mult. (361->99), div. (0->0), fcn. (450->14), ass. (0->51)
	t320 = sin(qJ(4));
	t325 = cos(qJ(4));
	t329 = pkin(12) + pkin(4);
	t355 = -qJ(5) * t325 + t329 * t320 + pkin(10);
	t315 = sin(pkin(7));
	t317 = cos(pkin(7));
	t321 = sin(qJ(3));
	t326 = cos(qJ(3));
	t330 = pkin(5) + pkin(11);
	t352 = qJ(5) * t320 + t329 * t325 + pkin(3);
	t335 = t321 * t352 - t330 * t326;
	t299 = -t355 * t315 + t335 * t317;
	t304 = t330 * t321 + t326 * t352 + pkin(2);
	t322 = sin(qJ(2));
	t327 = cos(qJ(2));
	t354 = t299 * t327 + t304 * t322;
	t353 = t335 * t315 + t355 * t317 + pkin(9);
	t350 = t315 * t321;
	t349 = t315 * t322;
	t348 = t315 * t326;
	t347 = t315 * t327;
	t346 = t320 * t321;
	t345 = t321 * t322;
	t344 = t321 * t327;
	t343 = t322 * t326;
	t342 = t326 * t327;
	t308 = t315 * t346 - t325 * t317;
	t316 = sin(pkin(6));
	t318 = cos(pkin(6));
	t307 = t315 * t325 + t317 * t346;
	t334 = t307 * t327 + t320 * t343;
	t297 = -t316 * t308 + t334 * t318;
	t333 = t317 * t342 - t345;
	t301 = -t316 * t348 + t333 * t318;
	t319 = sin(qJ(6));
	t324 = cos(qJ(6));
	t336 = t297 * t319 - t324 * t301;
	t332 = t317 * t344 + t343;
	t328 = cos(qJ(1));
	t323 = sin(qJ(1));
	t310 = t317 * t343 + t344;
	t309 = -t317 * t345 + t342;
	t306 = t316 * t317 + t318 * t347;
	t303 = -t322 * t307 + t320 * t342;
	t302 = t333 * t316 + t318 * t348;
	t300 = -t316 * t350 + t332 * t318;
	t298 = t303 * t319 + t324 * t310;
	t296 = t318 * t308 + t334 * t316;
	t295 = -t299 * t322 + t304 * t327 + pkin(1);
	t294 = -t353 * t316 + t354 * t318;
	t1 = [t298 * t328 - t336 * t323, (-t297 * t323 + t328 * t303) * t324 - (t301 * t323 + t328 * t310) * t319, (-t300 * t323 + t328 * t309) * t325 + t320 * (t306 * t323 + t328 * t349), -t294 * t323 + t295 * t328 + 0; t298 * t323 + t336 * t328, (t297 * t324 + t319 * t301) * t328 - t323 * (-t303 * t324 + t319 * t310), (t300 * t325 - t320 * t306) * t328 + t323 * (t309 * t325 + t320 * t349), t294 * t328 + t295 * t323 + 0; t296 * t319 - t324 * t302, t296 * t324 + t319 * t302, (t332 * t316 + t318 * t350) * t325 - t320 * (t316 * t347 - t318 * t317), t354 * t316 + t353 * t318 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end