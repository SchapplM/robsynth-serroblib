% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRR15 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:34
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRR15_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR15_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:34:32
	% EndTime: 2020-11-04 22:34:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:34:32
	% EndTime: 2020-11-04 22:34:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t154 = cos(qJ(1));
	t153 = sin(qJ(1));
	t1 = [t154, -t153, 0, 0; t153, t154, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:34:32
	% EndTime: 2020-11-04 22:34:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t155 = sin(pkin(6));
	t158 = sin(qJ(1));
	t166 = t158 * t155;
	t157 = sin(qJ(2));
	t165 = t158 * t157;
	t159 = cos(qJ(2));
	t164 = t158 * t159;
	t160 = cos(qJ(1));
	t163 = t160 * t155;
	t162 = t160 * t157;
	t161 = t160 * t159;
	t156 = cos(pkin(6));
	t1 = [-t156 * t165 + t161, -t156 * t164 - t162, t166, t160 * pkin(1) + pkin(9) * t166 + 0; t156 * t162 + t164, t156 * t161 - t165, -t163, t158 * pkin(1) - pkin(9) * t163 + 0; t155 * t157, t155 * t159, t156, t156 * pkin(9) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:34:32
	% EndTime: 2020-11-04 22:34:32
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->32), mult. (104->58), div. (0->0), fcn. (138->10), ass. (0->29)
	t173 = sin(pkin(7));
	t176 = cos(pkin(6));
	t195 = t173 * t176;
	t181 = cos(qJ(2));
	t194 = t173 * t181;
	t174 = sin(pkin(6));
	t175 = cos(pkin(7));
	t193 = t174 * t175;
	t192 = t176 * t181;
	t177 = sin(qJ(3));
	t178 = sin(qJ(2));
	t191 = t177 * t178;
	t190 = t177 * t181;
	t180 = cos(qJ(3));
	t189 = t178 * t180;
	t179 = sin(qJ(1));
	t188 = t179 * t178;
	t187 = t180 * t181;
	t182 = cos(qJ(1));
	t186 = t182 * t178;
	t185 = t182 * t181;
	t184 = -pkin(2) * t178 + pkin(10) * t194;
	t168 = -t173 * t174 + t175 * t192;
	t183 = t168 * t177 + t176 * t189;
	t172 = t175 * pkin(10) + pkin(9);
	t170 = t173 * t178 * pkin(10) + pkin(2) * t181 + pkin(1);
	t169 = -t175 * t191 + t187;
	t167 = t174 * t172 + t184 * t176;
	t1 = [t182 * t169 - t183 * t179, (-t168 * t179 - t175 * t186) * t180 - (-t176 * t188 + t185) * t177, (t179 * t192 + t186) * t173 + t179 * t193, t167 * t179 + t170 * t182 + 0; t179 * t169 + t183 * t182, (t168 * t180 - t176 * t191) * t182 - t179 * (t175 * t189 + t190), -(t176 * t185 - t188) * t173 - t182 * t193, -t167 * t182 + t170 * t179 + 0; t177 * t195 + (t175 * t190 + t189) * t174, t180 * t195 + (t175 * t187 - t191) * t174, -t174 * t194 + t176 * t175, t172 * t176 - t184 * t174 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:34:32
	% EndTime: 2020-11-04 22:34:33
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (70->37), mult. (149->63), div. (0->0), fcn. (183->10), ass. (0->32)
	t204 = sin(pkin(7));
	t206 = cos(pkin(7));
	t208 = sin(qJ(3));
	t211 = cos(qJ(3));
	t215 = pkin(3) * t208 - qJ(4) * t211;
	t198 = t204 * pkin(10) - t215 * t206;
	t201 = pkin(3) * t211 + qJ(4) * t208 + pkin(2);
	t209 = sin(qJ(2));
	t212 = cos(qJ(2));
	t231 = -t198 * t212 + t201 * t209;
	t230 = t206 * pkin(10) + t215 * t204 + pkin(9);
	t205 = sin(pkin(6));
	t227 = t204 * t205;
	t207 = cos(pkin(6));
	t226 = t204 * t207;
	t225 = t205 * t206;
	t224 = t207 * t212;
	t223 = t208 * t209;
	t222 = t208 * t212;
	t221 = t209 * t211;
	t210 = sin(qJ(1));
	t220 = t210 * t209;
	t219 = t211 * t212;
	t213 = cos(qJ(1));
	t218 = t213 * t209;
	t217 = t213 * t212;
	t199 = t206 * t224 - t227;
	t214 = t199 * t208 + t207 * t221;
	t200 = -t206 * t223 + t219;
	t197 = t198 * t209 + t201 * t212 + pkin(1);
	t196 = t205 * t230 - t231 * t207;
	t1 = [(t210 * t224 + t218) * t204 + t210 * t225, -t213 * t200 + t214 * t210, (t199 * t210 + t206 * t218) * t211 + (-t207 * t220 + t217) * t208, t196 * t210 + t197 * t213 + 0; -(t207 * t217 - t220) * t204 - t213 * t225, -t210 * t200 - t214 * t213, (-t199 * t211 + t207 * t223) * t213 + t210 * (t206 * t221 + t222), -t196 * t213 + t197 * t210 + 0; t207 * t206 - t212 * t227, -t208 * t226 + (-t206 * t222 - t221) * t205, -t211 * t226 + (-t206 * t219 + t223) * t205, t231 * t205 + t230 * t207 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:34:33
	% EndTime: 2020-11-04 22:34:33
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (112->45), mult. (210->75), div. (0->0), fcn. (265->12), ass. (0->40)
	t246 = sin(pkin(7));
	t248 = cos(pkin(7));
	t258 = pkin(4) + pkin(10);
	t251 = sin(qJ(3));
	t255 = cos(qJ(3));
	t259 = pkin(3) + pkin(11);
	t263 = qJ(4) * t255 - t251 * t259;
	t274 = t263 * t246 - t258 * t248 - pkin(9);
	t247 = sin(pkin(6));
	t273 = t246 * t247;
	t252 = sin(qJ(2));
	t272 = t246 * t252;
	t271 = t246 * t255;
	t249 = cos(pkin(6));
	t256 = cos(qJ(2));
	t270 = t249 * t256;
	t269 = t251 * t252;
	t268 = t251 * t256;
	t267 = t252 * t255;
	t266 = t255 * t256;
	t265 = t248 * t270;
	t234 = -t247 * t271 - t249 * t269 + t255 * t265;
	t238 = t246 * t270 + t247 * t248;
	t250 = sin(qJ(5));
	t254 = cos(qJ(5));
	t262 = t234 * t250 + t254 * t238;
	t237 = t246 * t258 + t263 * t248;
	t243 = qJ(4) * t251 + t259 * t255 + pkin(2);
	t261 = t237 * t256 - t243 * t252;
	t260 = (t265 - t273) * t251 + t249 * t267;
	t257 = cos(qJ(1));
	t253 = sin(qJ(1));
	t242 = t248 * t267 + t268;
	t241 = -t248 * t269 + t266;
	t239 = -t249 * t248 + t256 * t273;
	t236 = t242 * t250 + t254 * t272;
	t235 = -t249 * t271 + (-t248 * t266 + t269) * t247;
	t233 = t237 * t252 + t243 * t256 + pkin(1);
	t232 = -t247 * t274 + t261 * t249;
	t1 = [t236 * t257 + t262 * t253, (t234 * t253 + t257 * t242) * t254 - (t238 * t253 + t257 * t272) * t250, t257 * t241 - t260 * t253, t232 * t253 + t233 * t257 + 0; t236 * t253 - t262 * t257, (-t234 * t254 + t250 * t238) * t257 - t253 * (-t242 * t254 + t250 * t272), t253 * t241 + t260 * t257, -t232 * t257 + t233 * t253 + 0; t235 * t250 - t254 * t239, t235 * t254 + t250 * t239, t249 * t246 * t251 + (t248 * t268 + t267) * t247, -t261 * t247 - t274 * t249 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:34:33
	% EndTime: 2020-11-04 22:34:33
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (188->60), mult. (363->102), div. (0->0), fcn. (452->14), ass. (0->51)
	t298 = sin(qJ(5));
	t303 = cos(qJ(5));
	t330 = pkin(5) * t303 + pkin(12) * t298 + pkin(4) + pkin(10);
	t293 = sin(pkin(7));
	t295 = cos(pkin(7));
	t291 = -t298 * pkin(5) + pkin(12) * t303 - qJ(4);
	t299 = sin(qJ(3));
	t304 = cos(qJ(3));
	t307 = pkin(3) + pkin(11);
	t312 = t291 * t304 + t307 * t299;
	t280 = -t293 * t330 + t312 * t295;
	t300 = sin(qJ(2));
	t305 = cos(qJ(2));
	t311 = t291 * t299 - t307 * t304 - pkin(2);
	t329 = t280 * t305 - t311 * t300;
	t326 = t293 * t299;
	t325 = t293 * t300;
	t324 = t293 * t304;
	t323 = t295 * t303;
	t296 = cos(pkin(6));
	t322 = t296 * t305;
	t321 = t298 * t304;
	t320 = t299 * t300;
	t319 = t299 * t305;
	t318 = t300 * t304;
	t317 = t304 * t305;
	t294 = sin(pkin(6));
	t316 = t294 * t324;
	t315 = t298 * t320;
	t286 = t293 * t303 + t295 * t321;
	t278 = t286 * t322 + t294 * t323 - t296 * t315 - t298 * t316;
	t309 = t295 * t319 + t318;
	t282 = -t294 * t326 + t309 * t296;
	t297 = sin(qJ(6));
	t302 = cos(qJ(6));
	t313 = t278 * t297 + t302 * t282;
	t310 = t295 * t317 - t320;
	t308 = t312 * t293 + t330 * t295 + pkin(9);
	t306 = cos(qJ(1));
	t301 = sin(qJ(1));
	t288 = t295 * t318 + t319;
	t287 = t295 * t320 - t317;
	t285 = t293 * t322 + t294 * t295;
	t284 = t300 * t286 + t298 * t319;
	t283 = t309 * t294 + t296 * t326;
	t281 = t310 * t296 - t316;
	t279 = t284 * t297 + t302 * t287;
	t277 = (t293 * t321 - t323) * t296 + (t286 * t305 - t315) * t294;
	t276 = -t280 * t300 - t311 * t305 + pkin(1);
	t275 = -t308 * t294 + t329 * t296;
	t1 = [(t278 * t301 + t306 * t284) * t302 - (t282 * t301 + t306 * t287) * t297, -t306 * t279 - t313 * t301, (-t281 * t301 - t306 * t288) * t303 + (t285 * t301 + t306 * t325) * t298, -t275 * t301 + t276 * t306 + 0; (-t278 * t302 + t297 * t282) * t306 + (t284 * t302 - t297 * t287) * t301, -t279 * t301 + t313 * t306, (t281 * t303 - t298 * t285) * t306 + t301 * (-t288 * t303 + t298 * t325), t275 * t306 + t276 * t301 + 0; -t277 * t302 + t283 * t297, t277 * t297 + t283 * t302, (t310 * t294 + t296 * t324) * t303 - t298 * (t294 * t305 * t293 - t296 * t295), t329 * t294 + t308 * t296 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end