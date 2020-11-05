% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:10
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:10:04
	% EndTime: 2020-11-04 20:10:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:10:04
	% EndTime: 2020-11-04 20:10:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t182 = cos(pkin(11));
	t181 = sin(pkin(11));
	t1 = [t182, -t181, 0, 0; t181, t182, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:10:04
	% EndTime: 2020-11-04 20:10:04
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t183 = sin(pkin(11));
	t184 = sin(pkin(5));
	t192 = t183 * t184;
	t185 = cos(pkin(11));
	t191 = t185 * t184;
	t186 = cos(pkin(5));
	t187 = sin(qJ(2));
	t190 = t186 * t187;
	t188 = cos(qJ(2));
	t189 = t186 * t188;
	t1 = [-t183 * t190 + t185 * t188, -t183 * t189 - t185 * t187, t192, t185 * pkin(1) + pkin(7) * t192 + 0; t183 * t188 + t185 * t190, -t183 * t187 + t185 * t189, -t191, t183 * pkin(1) - pkin(7) * t191 + 0; t184 * t187, t184 * t188, t186, t186 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:10:04
	% EndTime: 2020-11-04 20:10:04
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->35), mult. (112->66), div. (0->0), fcn. (146->10), ass. (0->30)
	t199 = sin(pkin(6));
	t221 = pkin(8) * t199;
	t198 = sin(pkin(11));
	t220 = t198 * pkin(2);
	t201 = cos(pkin(11));
	t219 = t201 * pkin(2);
	t203 = cos(pkin(5));
	t218 = t199 * t203;
	t207 = cos(qJ(2));
	t217 = t199 * t207;
	t202 = cos(pkin(6));
	t197 = t202 * pkin(8) + pkin(7);
	t200 = sin(pkin(5));
	t216 = t200 * t197;
	t215 = t200 * t202;
	t205 = sin(qJ(2));
	t214 = t202 * t205;
	t213 = t202 * t207;
	t212 = t203 * t205;
	t211 = t203 * t207;
	t210 = t198 * t221;
	t209 = t201 * t221;
	t208 = -t199 * t200 + t202 * t211;
	t206 = cos(qJ(3));
	t204 = sin(qJ(3));
	t196 = t198 * t207 + t201 * t212;
	t195 = t198 * t212 - t201 * t207;
	t194 = -t198 * t214 + t208 * t201;
	t193 = -t208 * t198 - t201 * t214;
	t1 = [t193 * t204 - t206 * t195, t193 * t206 + t204 * t195, (t198 * t211 + t201 * t205) * t199 + t198 * t215, (t203 * t210 + t219) * t207 + (-t203 * t220 + t209) * t205 + t198 * t216 + t201 * pkin(1) + 0; t194 * t204 + t196 * t206, t194 * t206 - t196 * t204, -(-t198 * t205 + t201 * t211) * t199 - t201 * t215, (-t203 * t209 + t220) * t207 + (t203 * t219 + t210) * t205 - t201 * t216 + t198 * pkin(1) + 0; t204 * t218 + (t204 * t213 + t205 * t206) * t200, t206 * t218 + (-t204 * t205 + t206 * t213) * t200, -t200 * t217 + t203 * t202, t197 * t203 + qJ(1) + 0 + (pkin(2) * t205 - pkin(8) * t217) * t200; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:10:04
	% EndTime: 2020-11-04 20:10:04
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (91->59), mult. (247->109), div. (0->0), fcn. (302->12), ass. (0->42)
	t233 = sin(pkin(5));
	t232 = sin(pkin(6));
	t235 = cos(pkin(6));
	t238 = sin(qJ(3));
	t241 = cos(qJ(3));
	t250 = pkin(3) * t238 - pkin(9) * t241;
	t247 = t235 * pkin(8) + t250 * t232 + pkin(7);
	t268 = t247 * t233;
	t239 = sin(qJ(2));
	t242 = cos(qJ(2));
	t266 = t232 * pkin(8);
	t246 = -t250 * t235 + t266;
	t249 = pkin(3) * t241 + pkin(9) * t238 + pkin(2);
	t267 = t249 * t239 - t246 * t242;
	t231 = sin(pkin(11));
	t265 = t231 * t235;
	t264 = t232 * t233;
	t263 = t233 * t235;
	t234 = cos(pkin(11));
	t236 = cos(pkin(5));
	t262 = t234 * t236;
	t261 = t235 * t238;
	t260 = t235 * t239;
	t259 = t235 * t242;
	t258 = t236 * t235;
	t257 = t236 * t238;
	t256 = t236 * t239;
	t255 = t236 * t241;
	t254 = t236 * t242;
	t253 = t234 * t258;
	t252 = t235 * t257;
	t251 = t238 * t264;
	t248 = t235 * t254 - t264;
	t245 = -(t231 * t252 - t234 * t241) * t242 - (t231 * t255 + t234 * t261) * t239 + t231 * t251;
	t244 = -(t231 * t241 + t234 * t252) * t242 - (-t231 * t261 + t234 * t255) * t239 + t234 * t251;
	t240 = cos(qJ(4));
	t237 = sin(qJ(4));
	t229 = t242 * t264 - t258;
	t224 = t232 * t257 + (t238 * t259 + t241 * t239) * t233;
	t223 = t234 * t263 + (-t231 * t239 + t234 * t254) * t232;
	t222 = t232 * t239 * t234 + (t232 * t254 + t263) * t231;
	t1 = [t237 * t222 + t245 * t240, t240 * t222 - t245 * t237, (t248 * t231 + t234 * t260) * t241 - t238 * (t231 * t256 - t234 * t242), 0 + (t246 * t239 + t249 * t242 + pkin(1)) * t234 + (-t267 * t236 + t268) * t231; -t237 * t223 - t244 * t240, -t240 * t223 + t244 * t237, (t231 * t260 - t248 * t234) * t241 + (t231 * t242 + t234 * t256) * t238, ((t231 * pkin(3) - pkin(9) * t253) * t241 + (pkin(3) * t253 + t231 * pkin(9)) * t238 - t262 * t266 + t231 * pkin(2)) * t242 + ((pkin(3) * t262 + pkin(9) * t265) * t241 + (-pkin(3) * t265 + pkin(9) * t262) * t238 + pkin(2) * t262 + t231 * t266) * t239 + t231 * pkin(1) + 0 - t234 * t268; t224 * t240 - t237 * t229, -t224 * t237 - t240 * t229, -t232 * t255 + (t238 * t239 - t241 * t259) * t233, t267 * t233 + t247 * t236 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:10:04
	% EndTime: 2020-11-04 20:10:05
	% DurationCPUTime: 0.48s
	% Computational Cost: add. (167->81), mult. (469->151), div. (0->0), fcn. (558->14), ass. (0->62)
	t291 = cos(pkin(6));
	t298 = cos(qJ(4));
	t328 = t298 * t291;
	t294 = sin(qJ(4));
	t331 = t294 * t291;
	t358 = -pkin(4) * t331 - t291 * pkin(8) + pkin(10) * t328 - pkin(7);
	t357 = t298 * pkin(4) + t294 * pkin(10) + pkin(3);
	t288 = sin(pkin(6));
	t295 = sin(qJ(3));
	t299 = cos(qJ(3));
	t349 = t299 * pkin(9);
	t356 = (t295 * t357 - t349) * t291 - (t294 * pkin(4) - t298 * pkin(10) + pkin(8)) * t288;
	t287 = sin(pkin(11));
	t296 = sin(qJ(2));
	t290 = cos(pkin(11));
	t339 = t290 * t288;
	t289 = sin(pkin(5));
	t342 = t289 * t291;
	t275 = t287 * t342 + t296 * t339;
	t336 = t291 * t296;
	t346 = t288 * t289;
	t278 = -t287 * t346 + t290 * t336;
	t279 = -t288 * t294 + t295 * t328;
	t300 = cos(qJ(2));
	t292 = cos(pkin(5));
	t329 = t296 * t299;
	t321 = t292 * t329;
	t327 = t298 * t299;
	t348 = t287 * t292;
	t353 = (t278 * t295 + t287 * t321) * t298 - (-t279 * t348 + t290 * t327) * t300 - t294 * t275;
	t276 = t287 * t336 + t289 * t339;
	t282 = t290 * t342;
	t338 = t290 * t292;
	t347 = t287 * t296;
	t352 = (-t276 * t295 + t290 * t321) * t298 + (t279 * t338 + t287 * t327) * t300 - t294 * (-t288 * t347 + t282);
	t304 = pkin(9) * t295 + t299 * t357 + pkin(2);
	t351 = t358 * t289 + t346 * t349;
	t345 = t288 * t292;
	t344 = t288 * t295;
	t341 = t289 * t296;
	t340 = t289 * t300;
	t337 = t291 * t295;
	t335 = t291 * t300;
	t334 = t292 * t295;
	t333 = t292 * t299;
	t332 = t292 * t300;
	t330 = t295 * t296;
	t324 = t289 * t344;
	t323 = t291 * t334;
	t322 = t291 * t333;
	t320 = t292 * t330;
	t313 = t287 * t324;
	t312 = t290 * t324;
	t303 = t304 * t290;
	t301 = t356 * t287;
	t297 = cos(qJ(5));
	t293 = sin(qJ(5));
	t274 = t288 * t333 + (t299 * t335 - t330) * t289;
	t271 = -t292 * (t298 * t344 + t331) + (-t279 * t300 - t296 * t327) * t289;
	t270 = (t287 * t322 + t290 * t295) * t300 + t278 * t299 - t287 * t320;
	t269 = (-t287 * t295 + t290 * t322) * t300 - t276 * t299 - t290 * t320;
	t1 = [t293 * t270 - t353 * t297, t297 * t270 + t353 * t293, ((-t287 * t323 + t290 * t299) * t300 + (-t287 * t333 - t290 * t337) * t296 + t313) * t294 - t298 * (t288 * t287 * t332 + t275), (-t292 * t301 + t303) * t300 + (-t290 * t356 - t304 * t348) * t296 + t357 * t313 + t290 * pkin(1) + 0 - t351 * t287; -t293 * t269 + t352 * t297, -t297 * t269 - t352 * t293, ((t287 * t299 + t290 * t323) * t300 + (-t287 * t337 + t290 * t333) * t296 - t312) * t294 + t298 * (t282 + (t290 * t332 - t347) * t288), (t304 * t287 + t356 * t338) * t300 + (t292 * t303 - t301) * t296 - t357 * t312 + t287 * pkin(1) + 0 + t351 * t290; -t271 * t297 - t293 * t274, t271 * t293 - t297 * t274, (t288 * t334 + (t295 * t335 + t329) * t289) * t294 + t298 * (t288 * t340 - t292 * t291), t356 * t340 + (pkin(9) * t341 + t345 * t357) * t295 + (-pkin(9) * t345 + t341 * t357) * t299 + pkin(2) * t341 + qJ(1) + 0 - t358 * t292; 0, 0, 0, 1;];
	Tc_mdh = t1;
end