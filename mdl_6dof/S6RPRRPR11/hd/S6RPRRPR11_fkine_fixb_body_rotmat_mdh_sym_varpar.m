% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPR11 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:50
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRPR11_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR11_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:07
	% EndTime: 2020-11-04 21:50:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:07
	% EndTime: 2020-11-04 21:50:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t166 = cos(qJ(1));
	t165 = sin(qJ(1));
	t1 = [t166, -t165, 0, 0; t165, t166, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:07
	% EndTime: 2020-11-04 21:50:08
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t167 = sin(pkin(12));
	t171 = sin(qJ(1));
	t178 = t171 * t167;
	t168 = sin(pkin(6));
	t177 = t171 * t168;
	t169 = cos(pkin(12));
	t176 = t171 * t169;
	t172 = cos(qJ(1));
	t175 = t172 * t167;
	t174 = t172 * t168;
	t173 = t172 * t169;
	t170 = cos(pkin(6));
	t1 = [-t170 * t178 + t173, -t170 * t176 - t175, t177, pkin(1) * t172 + qJ(2) * t177 + 0; t170 * t175 + t176, t170 * t173 - t178, -t174, pkin(1) * t171 - qJ(2) * t174 + 0; t168 * t167, t168 * t169, t170, qJ(2) * t170 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:08
	% EndTime: 2020-11-04 21:50:08
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->32), mult. (101->58), div. (0->0), fcn. (135->10), ass. (0->31)
	t187 = sin(pkin(7));
	t210 = pkin(9) * t187;
	t186 = sin(pkin(12));
	t192 = sin(qJ(3));
	t209 = t186 * t192;
	t194 = cos(qJ(3));
	t208 = t186 * t194;
	t188 = sin(pkin(6));
	t207 = t188 * t187;
	t190 = cos(pkin(7));
	t206 = t188 * t190;
	t191 = cos(pkin(6));
	t205 = t191 * t190;
	t204 = t191 * t192;
	t203 = t191 * t194;
	t189 = cos(pkin(12));
	t202 = t192 * t189;
	t193 = sin(qJ(1));
	t201 = t193 * t186;
	t200 = t194 * t189;
	t195 = cos(qJ(1));
	t199 = t195 * t186;
	t198 = t195 * t189;
	t182 = -t186 * pkin(2) + t189 * t210;
	t184 = t190 * pkin(9) + qJ(2);
	t197 = t182 * t191 + t188 * t184;
	t179 = t189 * t205 - t207;
	t196 = t179 * t192 + t186 * t203;
	t181 = t189 * pkin(2) + t186 * t210 + pkin(1);
	t180 = t190 * t209 - t200;
	t1 = [-t195 * t180 - t196 * t193, (-t179 * t193 - t190 * t199) * t194 + t192 * (t191 * t201 - t198), (t193 * t191 * t189 + t199) * t187 + t193 * t206, t181 * t195 + t197 * t193 + 0; -t193 * t180 + t196 * t195, (t179 * t194 - t186 * t204) * t195 - t193 * (t190 * t208 + t202), -(t191 * t198 - t201) * t187 - t195 * t206, t181 * t193 - t197 * t195 + 0; t187 * t204 + (t190 * t202 + t208) * t188, t187 * t203 + (t190 * t200 - t209) * t188, -t189 * t207 + t205, -t182 * t188 + t184 * t191 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:08
	% EndTime: 2020-11-04 21:50:08
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (91->53), mult. (213->97), div. (0->0), fcn. (268->12), ass. (0->38)
	t230 = cos(pkin(12));
	t231 = cos(pkin(7));
	t232 = cos(pkin(6));
	t246 = t232 * t231;
	t228 = sin(pkin(7));
	t229 = sin(pkin(6));
	t249 = t229 * t228;
	t218 = t230 * t246 - t249;
	t227 = sin(pkin(12));
	t234 = sin(qJ(3));
	t237 = cos(qJ(3));
	t245 = t232 * t237;
	t212 = t218 * t234 + t227 * t245;
	t247 = t232 * t228;
	t216 = t229 * t231 + t230 * t247;
	t233 = sin(qJ(4));
	t236 = cos(qJ(4));
	t256 = t212 * t233 + t236 * t216;
	t220 = t230 * t228 * pkin(9) - t227 * pkin(2);
	t248 = t230 * t231;
	t221 = -t227 * pkin(3) + pkin(10) * t248;
	t222 = pkin(3) * t248 + t227 * pkin(10);
	t224 = t231 * pkin(9) + qJ(2);
	t255 = (pkin(3) * t249 - t222 * t232) * t234 - (pkin(10) * t249 - t221 * t232) * t237 + t220 * t232 + t229 * t224;
	t254 = t227 * t228;
	t253 = t227 * t231;
	t252 = t227 * t234;
	t251 = t227 * t237;
	t238 = cos(qJ(1));
	t250 = t227 * t238;
	t243 = t237 * t230;
	t239 = (t229 * t248 + t247) * t234 + t229 * t251;
	t235 = sin(qJ(1));
	t219 = -t231 * t252 + t243;
	t215 = t230 * t249 - t246;
	t214 = t219 * t233 - t236 * t254;
	t211 = (t230 * pkin(3) + pkin(10) * t253) * t237 + (-pkin(3) * t253 + t230 * pkin(10)) * t234 + pkin(9) * t254 + t230 * pkin(2) + pkin(1);
	t1 = [(-t212 * t235 + t238 * t219) * t236 + (t216 * t235 + t228 * t250) * t233, -t214 * t238 + t256 * t235, (t218 * t235 + t231 * t250) * t237 - t234 * (t235 * t232 * t227 - t238 * t230), t211 * t238 + t255 * t235 + 0; (t212 * t236 - t233 * t216) * t238 + t235 * (t219 * t236 + t233 * t254), -t214 * t235 - t256 * t238, (-t218 * t237 + t232 * t252) * t238 + t235 * (t234 * t230 + t231 * t251), t211 * t235 - t255 * t238 + 0; -t233 * t215 + t239 * t236, -t236 * t215 - t239 * t233, -t228 * t245 + (-t231 * t243 + t252) * t229, (-pkin(10) * t247 - t221 * t229) * t237 + (pkin(3) * t247 + t222 * t229) * t234 - t220 * t229 + t224 * t232 + 0 + pkin(8); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:08
	% EndTime: 2020-11-04 21:50:08
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (167->68), mult. (407->115), div. (0->0), fcn. (532->14), ass. (0->49)
	t287 = cos(pkin(12));
	t288 = cos(pkin(7));
	t289 = cos(pkin(6));
	t303 = t289 * t288;
	t284 = sin(pkin(7));
	t285 = sin(pkin(6));
	t306 = t285 * t284;
	t273 = t287 * t303 - t306;
	t283 = sin(pkin(12));
	t291 = sin(qJ(3));
	t294 = cos(qJ(3));
	t302 = t289 * t294;
	t267 = t273 * t291 + t283 * t302;
	t304 = t289 * t284;
	t271 = t285 * t288 + t287 * t304;
	t290 = sin(qJ(4));
	t293 = cos(qJ(4));
	t313 = t267 * t290 + t293 * t271;
	t275 = t287 * t284 * pkin(9) - t283 * pkin(2);
	t305 = t287 * t288;
	t276 = -t283 * pkin(3) + pkin(10) * t305;
	t277 = pkin(3) * t305 + t283 * pkin(10);
	t279 = t288 * pkin(9) + qJ(2);
	t312 = (pkin(3) * t306 - t277 * t289) * t291 - (pkin(10) * t306 - t276 * t289) * t294 + t275 * t289 + t285 * t279;
	t311 = t283 * t284;
	t310 = t283 * t288;
	t309 = t283 * t291;
	t308 = t283 * t294;
	t295 = cos(qJ(1));
	t307 = t283 * t295;
	t300 = t294 * t287;
	t296 = (t285 * t305 + t304) * t291 + t285 * t308;
	t292 = sin(qJ(1));
	t286 = cos(pkin(13));
	t282 = sin(pkin(13));
	t274 = -t288 * t309 + t300;
	t270 = t287 * t306 - t303;
	t269 = t274 * t290 - t293 * t311;
	t266 = -t284 * t302 + (-t288 * t300 + t309) * t285;
	t265 = (t287 * pkin(3) + pkin(10) * t310) * t294 + (-pkin(3) * t310 + t287 * pkin(10)) * t291 + pkin(9) * t311 + t287 * pkin(2) + pkin(1);
	t264 = (-t273 * t294 + t289 * t309) * t295 + t292 * (t291 * t287 + t288 * t308);
	t263 = (t273 * t292 + t288 * t307) * t294 - t291 * (t292 * t289 * t283 - t295 * t287);
	t262 = -t293 * t270 - t296 * t290;
	t261 = -t290 * t270 + t296 * t293;
	t260 = (-t267 * t292 + t295 * t274) * t293 + (t271 * t292 + t284 * t307) * t290;
	t259 = -t269 * t292 - t313 * t295;
	t258 = (t267 * t293 - t290 * t271) * t295 + t292 * (t274 * t293 + t290 * t311);
	t257 = -t269 * t295 + t313 * t292;
	t1 = [t260 * t286 + t263 * t282, -t260 * t282 + t263 * t286, -t257, t260 * pkin(4) - t257 * qJ(5) + t265 * t295 + t312 * t292 + 0; t258 * t286 + t264 * t282, -t258 * t282 + t264 * t286, -t259, t258 * pkin(4) - t259 * qJ(5) + t265 * t292 - t312 * t295 + 0; t261 * t286 + t266 * t282, -t261 * t282 + t266 * t286, -t262, t261 * pkin(4) - t262 * qJ(5) + (-pkin(10) * t304 - t276 * t285) * t294 + (pkin(3) * t304 + t277 * t285) * t291 - t275 * t285 + t279 * t289 + 0 + pkin(8); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:08
	% EndTime: 2020-11-04 21:50:08
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (222->73), mult. (445->118), div. (0->0), fcn. (546->16), ass. (0->59)
	t334 = sin(pkin(13)) * pkin(5) + pkin(10);
	t341 = sin(pkin(12));
	t344 = cos(pkin(12));
	t345 = cos(pkin(7));
	t372 = t344 * t345;
	t324 = -t341 * pkin(3) + t334 * t372;
	t325 = pkin(3) * t372 + t341 * t334;
	t346 = cos(pkin(6));
	t369 = t346 * t345;
	t342 = sin(pkin(7));
	t343 = sin(pkin(6));
	t374 = t343 * t342;
	t329 = t344 * t369 - t374;
	t331 = t344 * t342 * pkin(9) - t341 * pkin(2);
	t333 = t345 * pkin(9) + qJ(2);
	t335 = cos(pkin(13)) * pkin(5) + pkin(4);
	t347 = qJ(5) + pkin(11);
	t349 = sin(qJ(3));
	t352 = cos(qJ(3));
	t351 = cos(qJ(4));
	t348 = sin(qJ(4));
	t368 = t347 * t348;
	t360 = t335 * t351 + t368;
	t370 = t346 * t342;
	t327 = t343 * t345 + t344 * t370;
	t366 = t351 * t327;
	t367 = t348 * t327;
	t373 = t343 * t352;
	t376 = t341 * t351;
	t382 = ((t335 * t376 + t341 * t368 - t324) * t352 + t325 * t349 - t331) * t346 + t342 * t334 * t373 - t343 * t333 - t335 * t367 + t347 * t366 - (pkin(3) * t374 - t360 * t329) * t349;
	t375 = t341 * t352;
	t357 = t329 * t349 + t346 * t375;
	t381 = t357 * t348 + t366;
	t377 = t341 * t349;
	t371 = t345 * t349;
	t365 = t352 * t344;
	t364 = t341 * t373;
	t316 = t357 * t351 - t367;
	t322 = -t351 * t365 + t341 * (-t342 * t348 + t351 * t371);
	t350 = sin(qJ(1));
	t353 = cos(qJ(1));
	t362 = t316 * t353 - t350 * t322;
	t361 = t316 * t350 + t353 * t322;
	t359 = pkin(3) + t360;
	t328 = t343 * t372 + t370;
	t358 = t328 * t349 + t364;
	t340 = pkin(13) + qJ(6);
	t337 = cos(t340);
	t336 = sin(t340);
	t330 = t349 * t344 + t345 * t375;
	t326 = t344 * t374 - t369;
	t321 = (-t341 * t371 + t365) * t348 - t342 * t376;
	t320 = t329 * t352 - t346 * t377;
	t319 = t328 * t352 - t343 * t377;
	t318 = t320 * t353 - t350 * t330;
	t317 = t320 * t350 + t353 * t330;
	t315 = -t326 * t348 + t358 * t351;
	t314 = pkin(1) + (t334 * t349 + t359 * t352 + pkin(2)) * t344 + ((t335 * t348 - t347 * t351 + pkin(9)) * t342 + (t334 * t352 - t359 * t349) * t345) * t341;
	t1 = [t317 * t336 - t361 * t337, t337 * t317 + t361 * t336, t321 * t353 - t381 * t350, t314 * t353 - t382 * t350 + 0; -t318 * t336 + t362 * t337, -t337 * t318 - t362 * t336, t321 * t350 + t381 * t353, t314 * t350 + t382 * t353 + 0; t315 * t337 - t319 * t336, -t315 * t336 - t319 * t337, t351 * t326 + t358 * t348, (pkin(3) * t370 + t325 * t343 + t360 * t328) * t349 + (t326 * t347 + t335 * t364) * t351 + (-t335 * t326 + t347 * t364) * t348 + (-t324 * t343 - t334 * t370) * t352 - t331 * t343 + t333 * t346 + 0 + pkin(8); 0, 0, 0, 1;];
	Tc_mdh = t1;
end