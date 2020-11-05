% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:11
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:04
	% EndTime: 2020-11-04 21:11:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:04
	% EndTime: 2020-11-04 21:11:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t150 = cos(pkin(12));
	t149 = sin(pkin(12));
	t1 = [t150, -t149, 0, 0; t149, t150, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:04
	% EndTime: 2020-11-04 21:11:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t151 = sin(pkin(12));
	t152 = sin(pkin(6));
	t160 = t151 * t152;
	t153 = cos(pkin(12));
	t159 = t153 * t152;
	t154 = cos(pkin(6));
	t155 = sin(qJ(2));
	t158 = t154 * t155;
	t156 = cos(qJ(2));
	t157 = t154 * t156;
	t1 = [-t151 * t158 + t153 * t156, -t151 * t157 - t153 * t155, t160, t153 * pkin(1) + pkin(8) * t160 + 0; t151 * t156 + t153 * t158, -t151 * t155 + t153 * t157, -t159, t151 * pkin(1) - pkin(8) * t159 + 0; t152 * t155, t152 * t156, t154, t154 * pkin(8) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:04
	% EndTime: 2020-11-04 21:11:05
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->35), mult. (112->66), div. (0->0), fcn. (146->10), ass. (0->30)
	t167 = sin(pkin(7));
	t189 = pkin(9) * t167;
	t166 = sin(pkin(12));
	t188 = t166 * pkin(2);
	t169 = cos(pkin(12));
	t187 = t169 * pkin(2);
	t171 = cos(pkin(6));
	t186 = t167 * t171;
	t175 = cos(qJ(2));
	t185 = t167 * t175;
	t170 = cos(pkin(7));
	t165 = t170 * pkin(9) + pkin(8);
	t168 = sin(pkin(6));
	t184 = t168 * t165;
	t183 = t168 * t170;
	t173 = sin(qJ(2));
	t182 = t170 * t173;
	t181 = t170 * t175;
	t180 = t171 * t173;
	t179 = t171 * t175;
	t178 = t166 * t189;
	t177 = t169 * t189;
	t176 = -t167 * t168 + t170 * t179;
	t174 = cos(qJ(3));
	t172 = sin(qJ(3));
	t164 = t166 * t175 + t169 * t180;
	t163 = t166 * t180 - t169 * t175;
	t162 = -t166 * t182 + t176 * t169;
	t161 = -t176 * t166 - t169 * t182;
	t1 = [t161 * t172 - t174 * t163, t161 * t174 + t172 * t163, (t166 * t179 + t169 * t173) * t167 + t166 * t183, (t171 * t178 + t187) * t175 + (-t171 * t188 + t177) * t173 + t166 * t184 + t169 * pkin(1) + 0; t162 * t172 + t164 * t174, t162 * t174 - t164 * t172, -(-t166 * t173 + t169 * t179) * t167 - t169 * t183, (-t171 * t177 + t188) * t175 + (t171 * t187 + t178) * t173 - t169 * t184 + t166 * pkin(1) + 0; t172 * t186 + (t172 * t181 + t173 * t174) * t168, t174 * t186 + (-t172 * t173 + t174 * t181) * t168, -t168 * t185 + t171 * t170, t165 * t171 + qJ(1) + 0 + (pkin(2) * t173 - pkin(9) * t185) * t168; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:05
	% EndTime: 2020-11-04 21:11:05
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (125->43), mult. (130->65), div. (0->0), fcn. (168->16), ass. (0->41)
	t211 = qJ(3) + pkin(13);
	t207 = pkin(7) + t211;
	t200 = sin(t207);
	t208 = pkin(7) - t211;
	t201 = sin(t208);
	t196 = -t200 + t201;
	t233 = -t196 / 0.2e1;
	t202 = cos(t207);
	t203 = cos(t208);
	t198 = t202 + t203;
	t232 = t198 / 0.2e1;
	t217 = cos(pkin(6));
	t231 = t217 / 0.2e1;
	t230 = pkin(3) * sin(qJ(3));
	t213 = sin(pkin(7));
	t216 = cos(pkin(7));
	t218 = pkin(9) + qJ(4);
	t229 = t213 * t230 + t216 * t218 + pkin(8);
	t212 = sin(pkin(12));
	t214 = sin(pkin(6));
	t228 = t212 * t214;
	t215 = cos(pkin(12));
	t227 = t214 * t215;
	t226 = t214 * t216;
	t220 = sin(qJ(2));
	t225 = t217 * t220;
	t221 = cos(qJ(2));
	t224 = t217 * t221;
	t223 = t228 / 0.2e1;
	t222 = -t227 / 0.2e1;
	t210 = cos(t211);
	t209 = sin(t211);
	t206 = cos(qJ(3)) * pkin(3) + pkin(2);
	t199 = t203 - t202;
	t197 = t201 + t200;
	t195 = -t213 * t218 + t216 * t230;
	t193 = t212 * t221 + t215 * t225;
	t192 = -t212 * t220 + t215 * t224;
	t191 = t212 * t224 + t215 * t220;
	t190 = t212 * t225 - t215 * t221;
	t1 = [-t190 * t210 + t191 * t196 / 0.2e1 + t199 * t223, t209 * t190 - t191 * t198 / 0.2e1 + t197 * t223, t191 * t213 + t212 * t226, t215 * pkin(1) - t190 * t206 - t191 * t195 + t229 * t228 + 0; t192 * t233 + t193 * t210 + t199 * t222, t192 * t232 - t193 * t209 + t197 * t222, -t192 * t213 - t215 * t226, t212 * pkin(1) + t192 * t195 + t193 * t206 - t229 * t227 + 0; t199 * t231 + (t220 * t210 + t221 * t233) * t214, t197 * t231 + (-t220 * t209 + t221 * t232) * t214, -t214 * t221 * t213 + t217 * t216, qJ(1) + 0 + t229 * t217 + (t195 * t221 + t220 * t206) * t214; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:05
	% EndTime: 2020-11-04 21:11:05
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (183->65), mult. (301->108), div. (0->0), fcn. (356->20), ass. (0->55)
	t264 = sin(qJ(2));
	t267 = cos(qJ(2));
	t259 = cos(pkin(7));
	t255 = sin(pkin(7));
	t261 = pkin(9) + qJ(4);
	t288 = t255 * t261;
	t253 = sin(pkin(13));
	t257 = cos(pkin(13));
	t244 = t257 * pkin(4) + t253 * pkin(10) + pkin(3);
	t245 = -t253 * pkin(4) + t257 * pkin(10);
	t263 = sin(qJ(3));
	t266 = cos(qJ(3));
	t294 = -t244 * t263 + t245 * t266;
	t268 = t294 * t259 + t288;
	t273 = t244 * t266 + t245 * t263 + pkin(2);
	t297 = t273 * t264 - t268 * t267;
	t269 = -t255 * t294 + t259 * t261 + pkin(8);
	t252 = qJ(3) + pkin(13);
	t248 = pkin(7) + t252;
	t249 = pkin(7) - t252;
	t242 = sin(t249) + sin(t248);
	t293 = -t242 / 0.2e1;
	t243 = cos(t248) + cos(t249);
	t292 = -t243 / 0.2e1;
	t254 = sin(pkin(12));
	t289 = t254 * t259;
	t256 = sin(pkin(6));
	t258 = cos(pkin(12));
	t287 = t256 * t258;
	t286 = t256 * t259;
	t285 = t256 * t267;
	t260 = cos(pkin(6));
	t284 = t258 * t260;
	t283 = t258 * t264;
	t282 = t259 * t264;
	t281 = t260 * t255;
	t280 = t260 * t264;
	t279 = t260 * t267;
	t278 = t244 * t284;
	t277 = t245 * t284;
	t240 = t254 * t280 - t258 * t267;
	t250 = sin(t252);
	t251 = cos(t252);
	t270 = -t255 * t256 + t259 * t279;
	t276 = (t270 * t254 + t258 * t282) * t250 + t240 * t251;
	t241 = t254 * t267 + t258 * t280;
	t275 = (-t254 * t282 + t270 * t258) * t250 + t241 * t251;
	t272 = t251 * t256 * t264 + (t259 * t285 + t281) * t250;
	t271 = -t254 * t264 + t258 * t279;
	t265 = cos(qJ(5));
	t262 = sin(qJ(5));
	t238 = t255 * t285 - t260 * t259;
	t236 = t271 * t255 + t258 * t286;
	t235 = t255 * t283 + (t255 * t279 + t286) * t254;
	t1 = [t262 * t235 - t276 * t265, t265 * t235 + t276 * t262, -t250 * t240 + (t254 * t279 + t283) * t243 / 0.2e1 + t254 * t256 * t293, 0 + (t268 * t264 + t273 * t267 + pkin(1)) * t258 + (t269 * t256 - t297 * t260) * t254; -t262 * t236 + t275 * t265, -t265 * t236 - t275 * t262, t241 * t250 + t271 * t292 + t242 * t287 / 0.2e1, ((t254 * t245 + t259 * t278) * t263 + (t254 * t244 - t259 * t277) * t266 - t261 * t258 * t281 + t254 * pkin(2)) * t267 + ((-t244 * t289 + t277) * t263 + (t245 * t289 + t278) * t266 + pkin(2) * t284 + t254 * t288) * t264 + t254 * pkin(1) + 0 - t269 * t287; -t262 * t238 + t272 * t265, -t265 * t238 - t272 * t262, t260 * t293 + (t264 * t250 + t267 * t292) * t256, t297 * t256 + t269 * t260 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:11:05
	% EndTime: 2020-11-04 21:11:05
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (310->78), mult. (518->126), div. (0->0), fcn. (643->22), ass. (0->66)
	t338 = sin(qJ(2));
	t342 = cos(qJ(2));
	t332 = cos(pkin(7));
	t328 = sin(pkin(7));
	t334 = pkin(9) + qJ(4);
	t363 = t328 * t334;
	t326 = sin(pkin(13));
	t330 = cos(pkin(13));
	t317 = pkin(4) * t330 + pkin(10) * t326 + pkin(3);
	t318 = -pkin(4) * t326 + pkin(10) * t330;
	t337 = sin(qJ(3));
	t341 = cos(qJ(3));
	t369 = -t317 * t337 + t318 * t341;
	t343 = t332 * t369 + t363;
	t348 = t317 * t341 + t318 * t337 + pkin(2);
	t372 = t348 * t338 - t343 * t342;
	t344 = -t328 * t369 + t332 * t334 + pkin(8);
	t325 = qJ(3) + pkin(13);
	t321 = pkin(7) + t325;
	t322 = pkin(7) - t325;
	t315 = sin(t322) + sin(t321);
	t368 = -t315 / 0.2e1;
	t316 = cos(t321) + cos(t322);
	t367 = -t316 / 0.2e1;
	t327 = sin(pkin(12));
	t364 = t327 * t332;
	t329 = sin(pkin(6));
	t331 = cos(pkin(12));
	t362 = t329 * t331;
	t361 = t329 * t332;
	t360 = t329 * t342;
	t333 = cos(pkin(6));
	t359 = t331 * t333;
	t358 = t331 * t338;
	t357 = t332 * t338;
	t356 = t333 * t328;
	t355 = t333 * t338;
	t354 = t333 * t342;
	t353 = t317 * t359;
	t352 = t318 * t359;
	t313 = t327 * t355 - t331 * t342;
	t323 = sin(t325);
	t324 = cos(t325);
	t345 = -t328 * t329 + t332 * t354;
	t351 = (t327 * t345 + t331 * t357) * t323 + t313 * t324;
	t314 = t327 * t342 + t331 * t355;
	t350 = (-t327 * t357 + t331 * t345) * t323 + t314 * t324;
	t347 = t324 * t329 * t338 + (t332 * t360 + t356) * t323;
	t346 = -t327 * t338 + t331 * t354;
	t340 = cos(qJ(5));
	t339 = cos(qJ(6));
	t336 = sin(qJ(5));
	t335 = sin(qJ(6));
	t311 = t328 * t360 - t332 * t333;
	t309 = t328 * t346 + t331 * t361;
	t308 = t328 * t358 + (t328 * t354 + t361) * t327;
	t306 = t333 * t368 + (t338 * t323 + t342 * t367) * t329;
	t305 = -t336 * t311 + t340 * t347;
	t304 = t340 * t311 + t336 * t347;
	t303 = t314 * t323 + t346 * t367 + t315 * t362 / 0.2e1;
	t302 = -t323 * t313 + (t327 * t354 + t358) * t316 / 0.2e1 + t327 * t329 * t368;
	t301 = -t336 * t309 + t340 * t350;
	t300 = t336 * t308 - t340 * t351;
	t299 = t340 * t309 + t336 * t350;
	t298 = t340 * t308 + t336 * t351;
	t1 = [t300 * t339 + t302 * t335, -t300 * t335 + t302 * t339, -t298, t300 * pkin(5) - t298 * pkin(11) + 0 + (t338 * t343 + t342 * t348 + pkin(1)) * t331 + (t344 * t329 - t333 * t372) * t327; t301 * t339 + t303 * t335, -t301 * t335 + t303 * t339, t299, t301 * pkin(5) + t299 * pkin(11) + ((t318 * t327 + t332 * t353) * t337 + (t317 * t327 - t332 * t352) * t341 - t334 * t331 * t356 + t327 * pkin(2)) * t342 + ((-t317 * t364 + t352) * t337 + (t318 * t364 + t353) * t341 + pkin(2) * t359 + t327 * t363) * t338 + t327 * pkin(1) + 0 - t344 * t362; t305 * t339 + t306 * t335, -t305 * t335 + t306 * t339, t304, t305 * pkin(5) + t304 * pkin(11) + t329 * t372 + t344 * t333 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end