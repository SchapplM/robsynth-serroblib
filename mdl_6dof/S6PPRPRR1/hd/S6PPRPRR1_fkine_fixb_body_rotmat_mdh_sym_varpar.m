% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PPRPRR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:53
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PPRPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:53:15
	% EndTime: 2020-11-04 20:53:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:53:15
	% EndTime: 2020-11-04 20:53:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t149 = cos(pkin(11));
	t148 = sin(pkin(11));
	t1 = [t149, -t148, 0, 0; t148, t149, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:53:15
	% EndTime: 2020-11-04 20:53:16
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t151 = sin(pkin(11));
	t152 = sin(pkin(6));
	t159 = t151 * t152;
	t155 = cos(pkin(6));
	t158 = t151 * t155;
	t154 = cos(pkin(11));
	t157 = t154 * t152;
	t156 = t154 * t155;
	t153 = cos(pkin(12));
	t150 = sin(pkin(12));
	t1 = [-t150 * t158 + t154 * t153, -t154 * t150 - t153 * t158, t159, t154 * pkin(1) + qJ(2) * t159 + 0; t150 * t156 + t151 * t153, -t151 * t150 + t153 * t156, -t157, t151 * pkin(1) - qJ(2) * t157 + 0; t152 * t150, t152 * t153, t155, t155 * qJ(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:53:16
	% EndTime: 2020-11-04 20:53:16
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (44->30), mult. (101->54), div. (0->0), fcn. (135->10), ass. (0->29)
	t171 = sin(pkin(7));
	t188 = pkin(8) * t171;
	t169 = sin(pkin(12));
	t174 = cos(pkin(11));
	t187 = t169 * t174;
	t170 = sin(pkin(11));
	t186 = t170 * t169;
	t176 = cos(pkin(6));
	t185 = t171 * t176;
	t172 = sin(pkin(6));
	t184 = t172 * t171;
	t173 = cos(pkin(12));
	t175 = cos(pkin(7));
	t183 = t173 * t175;
	t182 = t174 * t176;
	t181 = t175 * t172;
	t180 = t176 * t175;
	t166 = -t169 * pkin(2) + t173 * t188;
	t167 = t175 * pkin(8) + qJ(2);
	t179 = t166 * t176 + t172 * t167;
	t178 = cos(qJ(3));
	t177 = sin(qJ(3));
	t165 = t173 * pkin(2) + t169 * t188 + pkin(1);
	t164 = t174 * t173 - t176 * t186;
	t163 = t169 * t182 + t170 * t173;
	t162 = t173 * t180 - t184;
	t161 = -t162 * t170 - t175 * t187;
	t160 = t162 * t174 - t175 * t186;
	t1 = [t161 * t177 + t164 * t178, t161 * t178 - t164 * t177, (t173 * t185 + t181) * t170 + t171 * t187, t165 * t174 + t179 * t170 + 0; t160 * t177 + t163 * t178, t160 * t178 - t163 * t177, (-t173 * t182 + t186) * t171 - t174 * t181, t165 * t170 - t179 * t174 + 0; t177 * t185 + (t169 * t178 + t177 * t183) * t172, t178 * t185 + (-t169 * t177 + t178 * t183) * t172, -t173 * t184 + t180, -t166 * t172 + t167 * t176 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:53:16
	% EndTime: 2020-11-04 20:53:16
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (125->47), mult. (131->72), div. (0->0), fcn. (170->16), ass. (0->47)
	t218 = cos(pkin(7));
	t220 = pkin(8) + qJ(4);
	t238 = t218 * t220 + qJ(2);
	t211 = qJ(3) + pkin(13);
	t207 = pkin(7) + t211;
	t201 = sin(t207);
	t208 = pkin(7) - t211;
	t202 = sin(t208);
	t195 = -t201 + t202;
	t237 = t195 / 0.2e1;
	t203 = cos(t207);
	t204 = cos(t208);
	t197 = t203 + t204;
	t236 = -t197 / 0.2e1;
	t219 = cos(pkin(6));
	t235 = t219 / 0.2e1;
	t221 = sin(qJ(3));
	t234 = pkin(3) * t221;
	t212 = sin(pkin(12));
	t213 = sin(pkin(11));
	t233 = t212 * t213;
	t217 = cos(pkin(11));
	t232 = t212 * t217;
	t215 = sin(pkin(6));
	t231 = t213 * t215;
	t214 = sin(pkin(7));
	t230 = t214 * t220;
	t229 = t215 * t217;
	t228 = t215 * t218;
	t216 = cos(pkin(12));
	t227 = t216 * t219;
	t226 = t217 * t219;
	t225 = t214 * t234 + t238;
	t224 = t231 / 0.2e1;
	t223 = -t229 / 0.2e1;
	t222 = cos(qJ(3));
	t210 = cos(t211);
	t209 = sin(t211);
	t206 = pkin(3) * t222 + pkin(2);
	t198 = t204 - t203;
	t196 = t202 + t201;
	t194 = t218 * t234 - t230;
	t192 = t212 * t226 + t213 * t216;
	t191 = -t216 * t226 + t233;
	t190 = t213 * t227 + t232;
	t189 = -t217 * t216 + t219 * t233;
	t1 = [-t189 * t210 + t190 * t237 + t198 * t224, t189 * t209 + t190 * t236 + t196 * t224, (t214 * t227 + t228) * t213 + t214 * t232, t217 * pkin(1) - t189 * t206 - t190 * t194 + t225 * t231 + 0; t191 * t237 + t192 * t210 + t198 * t223, t191 * t236 - t192 * t209 + t196 * t223, t191 * t214 - t217 * t228, t213 * pkin(1) - t191 * t194 + t192 * t206 - t225 * t229 + 0; t198 * t235 + (t212 * t210 - t216 * t195 / 0.2e1) * t215, t196 * t235 + (-t212 * t209 + t216 * t197 / 0.2e1) * t215, -t214 * t215 * t216 + t218 * t219, (pkin(2) * t212 - t216 * t230) * t215 + t238 * t219 + 0 + qJ(1) + ((t214 * t219 + t216 * t228) * t221 + t222 * t212 * t215) * pkin(3); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:53:16
	% EndTime: 2020-11-04 20:53:16
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (183->66), mult. (277->111), div. (0->0), fcn. (332->20), ass. (0->58)
	t267 = qJ(3) + pkin(13);
	t262 = pkin(7) + t267;
	t263 = pkin(7) - t267;
	t257 = cos(t262) + cos(t263);
	t301 = -t257 / 0.2e1;
	t269 = sin(pkin(12));
	t276 = cos(pkin(7));
	t300 = t269 * t276;
	t270 = sin(pkin(11));
	t299 = t270 * t269;
	t277 = cos(pkin(6));
	t298 = t270 * t277;
	t271 = sin(pkin(7));
	t278 = pkin(8) + qJ(4);
	t297 = t271 * t278;
	t256 = sin(t263) + sin(t262);
	t272 = sin(pkin(6));
	t296 = t272 * t256;
	t295 = t272 * t271;
	t274 = cos(pkin(12));
	t294 = t274 * t276;
	t275 = cos(pkin(11));
	t293 = t275 * t269;
	t292 = t275 * t277;
	t291 = t276 * t272;
	t290 = t277 * t271;
	t289 = t277 * t276;
	t247 = t269 * t298 - t275 * t274;
	t250 = t274 * t289 - t295;
	t264 = sin(t267);
	t265 = cos(t267);
	t288 = (t250 * t270 + t276 * t293) * t264 + t247 * t265;
	t252 = t269 * t292 + t270 * t274;
	t287 = (t250 * t275 - t276 * t299) * t264 + t252 * t265;
	t254 = -t269 * pkin(2) + t274 * t297;
	t260 = t276 * t278 + qJ(2);
	t286 = t254 * t277 + t272 * t260;
	t268 = sin(pkin(13));
	t273 = cos(pkin(13));
	t258 = t273 * pkin(4) + t268 * pkin(9) + pkin(3);
	t259 = -t268 * pkin(4) + t273 * pkin(9);
	t239 = t258 * t294 + t269 * t259;
	t285 = -t239 * t277 + t258 * t295;
	t240 = -t269 * t258 + t259 * t294;
	t284 = -t240 * t277 + t259 * t295;
	t283 = t265 * t269 * t272 + (t274 * t291 + t290) * t264;
	t282 = cos(qJ(3));
	t281 = cos(qJ(5));
	t280 = sin(qJ(3));
	t279 = sin(qJ(5));
	t253 = t274 * pkin(2) + t269 * t297 + pkin(1);
	t249 = t274 * t295 - t289;
	t248 = t274 * t290 + t291;
	t244 = t248 * t275 - t271 * t299;
	t243 = t248 * t270 + t271 * t293;
	t242 = t258 * t300 - t259 * t274;
	t241 = t258 * t274 + t259 * t300;
	t1 = [t243 * t279 - t288 * t281, t243 * t281 + t288 * t279, -t247 * t264 + (t274 * t298 + t293) * t257 / 0.2e1 - t270 * t296 / 0.2e1, (-t242 * t275 + t285 * t270) * t280 + (t275 * t241 - t284 * t270) * t282 + t286 * t270 + t253 * t275 + 0; -t279 * t244 + t287 * t281, -t244 * t281 - t287 * t279, t252 * t264 + (t274 * t292 - t299) * t301 + t275 * t296 / 0.2e1, (-t270 * t242 - t285 * t275) * t280 + (t270 * t241 + t284 * t275) * t282 - t286 * t275 + t253 * t270 + 0; -t279 * t249 + t283 * t281, -t281 * t249 - t283 * t279, -t277 * t256 / 0.2e1 + (t269 * t264 + t274 * t301) * t272, (t239 * t272 + t258 * t290) * t280 + (-t240 * t272 - t259 * t290) * t282 - t254 * t272 + t260 * t277 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:53:16
	% EndTime: 2020-11-04 20:53:16
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (310->79), mult. (482->129), div. (0->0), fcn. (607->22), ass. (0->69)
	t339 = qJ(3) + pkin(13);
	t334 = pkin(7) + t339;
	t335 = pkin(7) - t339;
	t329 = cos(t334) + cos(t335);
	t375 = -t329 / 0.2e1;
	t341 = sin(pkin(12));
	t348 = cos(pkin(7));
	t374 = t341 * t348;
	t342 = sin(pkin(11));
	t373 = t342 * t341;
	t349 = cos(pkin(6));
	t372 = t342 * t349;
	t343 = sin(pkin(7));
	t350 = pkin(8) + qJ(4);
	t371 = t343 * t350;
	t328 = sin(t335) + sin(t334);
	t344 = sin(pkin(6));
	t370 = t344 * t328;
	t369 = t344 * t343;
	t346 = cos(pkin(12));
	t368 = t346 * t348;
	t347 = cos(pkin(11));
	t367 = t347 * t341;
	t366 = t347 * t349;
	t365 = t348 * t344;
	t364 = t349 * t343;
	t363 = t349 * t348;
	t319 = t341 * t372 - t347 * t346;
	t322 = t346 * t363 - t369;
	t336 = sin(t339);
	t337 = cos(t339);
	t362 = (t322 * t342 + t348 * t367) * t336 + t319 * t337;
	t324 = t341 * t366 + t342 * t346;
	t361 = (t322 * t347 - t348 * t373) * t336 + t324 * t337;
	t326 = -t341 * pkin(2) + t346 * t371;
	t332 = t348 * t350 + qJ(2);
	t360 = t326 * t349 + t344 * t332;
	t340 = sin(pkin(13));
	t345 = cos(pkin(13));
	t330 = t345 * pkin(4) + t340 * pkin(9) + pkin(3);
	t331 = -t340 * pkin(4) + t345 * pkin(9);
	t311 = t330 * t368 + t341 * t331;
	t359 = -t311 * t349 + t330 * t369;
	t312 = -t341 * t330 + t331 * t368;
	t358 = -t312 * t349 + t331 * t369;
	t357 = t337 * t341 * t344 + (t346 * t365 + t364) * t336;
	t356 = cos(qJ(3));
	t355 = cos(qJ(5));
	t354 = cos(qJ(6));
	t353 = sin(qJ(3));
	t352 = sin(qJ(5));
	t351 = sin(qJ(6));
	t325 = t346 * pkin(2) + t341 * t371 + pkin(1);
	t321 = t346 * t369 - t363;
	t320 = t346 * t364 + t365;
	t316 = t320 * t347 - t343 * t373;
	t315 = t320 * t342 + t343 * t367;
	t314 = t330 * t374 - t331 * t346;
	t313 = t330 * t346 + t331 * t374;
	t310 = -t349 * t328 / 0.2e1 + (t341 * t336 + t346 * t375) * t344;
	t309 = -t352 * t321 + t357 * t355;
	t308 = t355 * t321 + t357 * t352;
	t307 = t324 * t336 + (t346 * t366 - t373) * t375 + t347 * t370 / 0.2e1;
	t306 = -t319 * t336 + (t346 * t372 + t367) * t329 / 0.2e1 - t342 * t370 / 0.2e1;
	t305 = -t352 * t316 + t361 * t355;
	t304 = t315 * t352 - t362 * t355;
	t303 = t316 * t355 + t361 * t352;
	t302 = t315 * t355 + t362 * t352;
	t1 = [t304 * t354 + t306 * t351, -t304 * t351 + t306 * t354, -t302, t304 * pkin(5) - t302 * pkin(10) + (-t314 * t347 + t359 * t342) * t353 + (t347 * t313 - t358 * t342) * t356 + t360 * t342 + t325 * t347 + 0; t305 * t354 + t307 * t351, -t305 * t351 + t307 * t354, t303, t305 * pkin(5) + t303 * pkin(10) + (-t342 * t314 - t359 * t347) * t353 + (t342 * t313 + t358 * t347) * t356 - t360 * t347 + t325 * t342 + 0; t309 * t354 + t310 * t351, -t309 * t351 + t310 * t354, t308, t309 * pkin(5) + t308 * pkin(10) + (t311 * t344 + t330 * t364) * t353 + (-t312 * t344 - t331 * t364) * t356 - t326 * t344 + t332 * t349 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end