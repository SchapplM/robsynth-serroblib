% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRR9 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:41
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:41:55
	% EndTime: 2020-11-04 21:41:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:41:55
	% EndTime: 2020-11-04 21:41:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t156 = cos(qJ(1));
	t155 = sin(qJ(1));
	t1 = [t156, -t155, 0, 0; t155, t156, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:41:55
	% EndTime: 2020-11-04 21:41:56
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t157 = sin(pkin(12));
	t161 = sin(qJ(1));
	t168 = t161 * t157;
	t158 = sin(pkin(6));
	t167 = t161 * t158;
	t159 = cos(pkin(12));
	t166 = t161 * t159;
	t162 = cos(qJ(1));
	t165 = t162 * t157;
	t164 = t162 * t158;
	t163 = t162 * t159;
	t160 = cos(pkin(6));
	t1 = [-t160 * t168 + t163, -t160 * t166 - t165, t167, t162 * pkin(1) + qJ(2) * t167 + 0; t160 * t165 + t166, t160 * t163 - t168, -t164, t161 * pkin(1) - qJ(2) * t164 + 0; t158 * t157, t158 * t159, t160, t160 * qJ(2) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:41:56
	% EndTime: 2020-11-04 21:41:56
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (44->32), mult. (101->58), div. (0->0), fcn. (135->10), ass. (0->31)
	t177 = sin(pkin(7));
	t200 = pkin(9) * t177;
	t176 = sin(pkin(12));
	t182 = sin(qJ(3));
	t199 = t176 * t182;
	t184 = cos(qJ(3));
	t198 = t176 * t184;
	t178 = sin(pkin(6));
	t197 = t178 * t177;
	t180 = cos(pkin(7));
	t196 = t178 * t180;
	t181 = cos(pkin(6));
	t195 = t181 * t180;
	t194 = t181 * t182;
	t193 = t181 * t184;
	t179 = cos(pkin(12));
	t192 = t182 * t179;
	t183 = sin(qJ(1));
	t191 = t183 * t176;
	t190 = t184 * t179;
	t185 = cos(qJ(1));
	t189 = t185 * t176;
	t188 = t185 * t179;
	t172 = -t176 * pkin(2) + t179 * t200;
	t174 = t180 * pkin(9) + qJ(2);
	t187 = t172 * t181 + t178 * t174;
	t169 = t179 * t195 - t197;
	t186 = t169 * t182 + t176 * t193;
	t171 = t179 * pkin(2) + t176 * t200 + pkin(1);
	t170 = t180 * t199 - t190;
	t1 = [-t185 * t170 - t186 * t183, (-t169 * t183 - t180 * t189) * t184 + t182 * (t181 * t191 - t188), (t183 * t181 * t179 + t189) * t177 + t183 * t196, t171 * t185 + t187 * t183 + 0; -t183 * t170 + t186 * t185, (t169 * t184 - t176 * t194) * t185 - t183 * (t180 * t198 + t192), -(t181 * t188 - t191) * t177 - t185 * t196, t171 * t183 - t187 * t185 + 0; t177 * t194 + (t180 * t192 + t198) * t178, t177 * t193 + (t180 * t190 - t199) * t178, -t179 * t197 + t195, -t172 * t178 + t174 * t181 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:41:56
	% EndTime: 2020-11-04 21:41:56
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (125->43), mult. (142->61), div. (0->0), fcn. (168->16), ass. (0->37)
	t238 = pkin(3) * sin(qJ(3));
	t223 = sin(pkin(6));
	t229 = sin(qJ(1));
	t237 = t223 * t229;
	t230 = cos(qJ(1));
	t236 = t223 * t230;
	t221 = sin(pkin(12));
	t235 = t229 * t221;
	t224 = cos(pkin(12));
	t234 = t229 * t224;
	t233 = t230 * t221;
	t232 = t230 * t224;
	t222 = sin(pkin(7));
	t225 = cos(pkin(7));
	t227 = pkin(9) + qJ(4);
	t231 = t222 * t238 + t225 * t227 + qJ(2);
	t220 = qJ(3) + pkin(13);
	t226 = cos(pkin(6));
	t219 = cos(t220);
	t218 = sin(t220);
	t217 = pkin(7) - t220;
	t216 = pkin(7) + t220;
	t215 = cos(qJ(3)) * pkin(3) + pkin(2);
	t214 = cos(t216);
	t213 = sin(t217);
	t212 = cos(t217) / 0.2e1;
	t211 = sin(t216) / 0.2e1;
	t210 = -t222 * t227 + t225 * t238;
	t208 = t212 - t214 / 0.2e1;
	t207 = t214 / 0.2e1 + t212;
	t206 = t213 / 0.2e1 + t211;
	t205 = t211 - t213 / 0.2e1;
	t204 = -t226 * t235 + t232;
	t203 = -t226 * t234 - t233;
	t202 = t226 * t233 + t234;
	t201 = t226 * t232 - t235;
	t1 = [t203 * t205 + t204 * t219 + t208 * t237, t203 * t207 - t204 * t218 + t206 * t237, -t203 * t222 + t225 * t237, t230 * pkin(1) + t203 * t210 + t204 * t215 + t231 * t237 + 0; t201 * t205 + t202 * t219 - t208 * t236, t201 * t207 - t202 * t218 - t206 * t236, -t201 * t222 - t225 * t236, t229 * pkin(1) + t201 * t210 + t202 * t215 - t231 * t236 + 0; t226 * t208 + (t205 * t224 + t219 * t221) * t223, t226 * t206 + (t207 * t224 - t218 * t221) * t223, -t223 * t224 * t222 + t226 * t225, pkin(8) + 0 + t231 * t226 + (t210 * t224 + t215 * t221) * t223; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:41:56
	% EndTime: 2020-11-04 21:41:56
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (183->62), mult. (273->101), div. (0->0), fcn. (328->20), ass. (0->54)
	t265 = sin(pkin(13));
	t269 = cos(pkin(13));
	t256 = t269 * pkin(4) + t265 * pkin(10) + pkin(3);
	t257 = -t265 * pkin(4) + t269 * pkin(10);
	t266 = sin(pkin(12));
	t270 = cos(pkin(12));
	t271 = cos(pkin(7));
	t293 = t270 * t271;
	t240 = t256 * t293 + t266 * t257;
	t241 = -t266 * t256 + t257 * t293;
	t267 = sin(pkin(7));
	t273 = pkin(9) + qJ(4);
	t296 = t267 * t273;
	t252 = -t266 * pkin(2) + t270 * t296;
	t258 = t271 * t273 + qJ(2);
	t268 = sin(pkin(6));
	t272 = cos(pkin(6));
	t275 = sin(qJ(3));
	t278 = cos(qJ(3));
	t294 = t268 * t267;
	t299 = (-t240 * t272 + t256 * t294) * t275 - (-t241 * t272 + t257 * t294) * t278 + t252 * t272 + t268 * t258;
	t264 = qJ(3) + pkin(13);
	t259 = pkin(7) + t264;
	t260 = pkin(7) - t264;
	t255 = cos(t259) + cos(t260);
	t298 = -t255 / 0.2e1;
	t297 = t266 * t271;
	t254 = sin(t260) + sin(t259);
	t295 = t268 * t254;
	t292 = t271 * t268;
	t291 = t272 * t267;
	t290 = t272 * t271;
	t276 = sin(qJ(1));
	t289 = t276 * t266;
	t288 = t276 * t270;
	t279 = cos(qJ(1));
	t287 = t279 * t266;
	t286 = t279 * t270;
	t249 = t270 * t290 - t294;
	t251 = -t272 * t289 + t286;
	t261 = sin(t264);
	t262 = cos(t264);
	t285 = (t249 * t276 + t271 * t287) * t261 - t251 * t262;
	t250 = t272 * t287 + t288;
	t284 = (t249 * t279 - t271 * t289) * t261 + t250 * t262;
	t280 = t262 * t266 * t268 + (t270 * t292 + t291) * t261;
	t277 = cos(qJ(5));
	t274 = sin(qJ(5));
	t247 = t270 * t291 + t292;
	t246 = t270 * t294 - t290;
	t243 = t247 * t279 - t267 * t289;
	t242 = t247 * t276 + t267 * t287;
	t239 = (-t256 * t297 + t257 * t270) * t275 + (t256 * t270 + t257 * t297) * t278 + t270 * pkin(2) + t266 * t296 + pkin(1);
	t1 = [t274 * t242 - t285 * t277, t277 * t242 + t285 * t274, t251 * t261 + (t272 * t288 + t287) * t255 / 0.2e1 - t276 * t295 / 0.2e1, t239 * t279 + t299 * t276 + 0; -t274 * t243 + t284 * t277, -t277 * t243 - t284 * t274, t250 * t261 + (t272 * t286 - t289) * t298 + t279 * t295 / 0.2e1, t239 * t276 - t299 * t279 + 0; -t274 * t246 + t280 * t277, -t277 * t246 - t280 * t274, -t272 * t254 / 0.2e1 + (t266 * t261 + t270 * t298) * t268, (t240 * t268 + t256 * t291) * t275 + (-t241 * t268 - t257 * t291) * t278 - t252 * t268 + t258 * t272 + 0 + pkin(8); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:41:56
	% EndTime: 2020-11-04 21:41:56
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (310->75), mult. (478->119), div. (0->0), fcn. (603->22), ass. (0->65)
	t335 = sin(pkin(13));
	t339 = cos(pkin(13));
	t326 = t339 * pkin(4) + t335 * pkin(10) + pkin(3);
	t327 = -t335 * pkin(4) + t339 * pkin(10);
	t336 = sin(pkin(12));
	t340 = cos(pkin(12));
	t341 = cos(pkin(7));
	t365 = t340 * t341;
	t310 = t326 * t365 + t336 * t327;
	t311 = -t336 * t326 + t327 * t365;
	t337 = sin(pkin(7));
	t343 = pkin(9) + qJ(4);
	t368 = t337 * t343;
	t322 = -t336 * pkin(2) + t340 * t368;
	t328 = t341 * t343 + qJ(2);
	t338 = sin(pkin(6));
	t342 = cos(pkin(6));
	t346 = sin(qJ(3));
	t350 = cos(qJ(3));
	t366 = t338 * t337;
	t371 = (-t310 * t342 + t326 * t366) * t346 - (-t311 * t342 + t327 * t366) * t350 + t322 * t342 + t338 * t328;
	t334 = qJ(3) + pkin(13);
	t329 = pkin(7) + t334;
	t330 = pkin(7) - t334;
	t325 = cos(t329) + cos(t330);
	t370 = -t325 / 0.2e1;
	t369 = t336 * t341;
	t324 = sin(t330) + sin(t329);
	t367 = t338 * t324;
	t364 = t341 * t338;
	t363 = t342 * t337;
	t362 = t342 * t341;
	t347 = sin(qJ(1));
	t361 = t347 * t336;
	t360 = t347 * t340;
	t351 = cos(qJ(1));
	t359 = t351 * t336;
	t358 = t351 * t340;
	t319 = t340 * t362 - t366;
	t321 = -t342 * t361 + t358;
	t331 = sin(t334);
	t332 = cos(t334);
	t357 = (t319 * t347 + t341 * t359) * t331 - t321 * t332;
	t320 = t342 * t359 + t360;
	t356 = (t319 * t351 - t341 * t361) * t331 + t320 * t332;
	t352 = t332 * t336 * t338 + (t340 * t364 + t363) * t331;
	t349 = cos(qJ(5));
	t348 = cos(qJ(6));
	t345 = sin(qJ(5));
	t344 = sin(qJ(6));
	t317 = t340 * t363 + t364;
	t316 = t340 * t366 - t362;
	t313 = t317 * t351 - t337 * t361;
	t312 = t317 * t347 + t337 * t359;
	t309 = -t342 * t324 / 0.2e1 + (t336 * t331 + t340 * t370) * t338;
	t308 = -t345 * t316 + t352 * t349;
	t307 = t349 * t316 + t352 * t345;
	t306 = t321 * t331 + (t342 * t360 + t359) * t325 / 0.2e1 - t347 * t367 / 0.2e1;
	t305 = t320 * t331 + (t342 * t358 - t361) * t370 + t351 * t367 / 0.2e1;
	t304 = (-t326 * t369 + t327 * t340) * t346 + (t326 * t340 + t327 * t369) * t350 + t340 * pkin(2) + t336 * t368 + pkin(1);
	t303 = -t345 * t313 + t356 * t349;
	t302 = t345 * t312 - t357 * t349;
	t301 = t349 * t313 + t356 * t345;
	t300 = t349 * t312 + t357 * t345;
	t1 = [t302 * t348 + t306 * t344, -t302 * t344 + t306 * t348, -t300, t302 * pkin(5) - t300 * pkin(11) + t304 * t351 + t371 * t347 + 0; t303 * t348 + t305 * t344, -t303 * t344 + t305 * t348, t301, t303 * pkin(5) + t301 * pkin(11) + t304 * t347 - t371 * t351 + 0; t308 * t348 + t309 * t344, -t308 * t344 + t309 * t348, t307, t308 * pkin(5) + t307 * pkin(11) + (t310 * t338 + t326 * t363) * t346 + (-t311 * t338 - t327 * t363) * t350 - t322 * t338 + t328 * t342 + 0 + pkin(8); 0, 0, 0, 1;];
	Tc_mdh = t1;
end