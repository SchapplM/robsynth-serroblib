% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6PRPRRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:44
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRP2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:35
	% EndTime: 2019-10-09 21:44:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:35
	% EndTime: 2019-10-09 21:44:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:35
	% EndTime: 2019-10-09 21:44:35
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(10));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(10)) * cos(pkin(6));
	t37 = t39 * t42 - t41 * t45;
	t34 = 0.1e1 / t37 ^ 2;
	t49 = qJD(2) * t34;
	t36 = t39 * t41 + t42 * t45;
	t33 = t36 ^ 2;
	t30 = t33 * t34 + 0.1e1;
	t46 = t37 * t49;
	t47 = t36 / t37 * t49;
	t48 = (t33 * t47 + t36 * t46) / t30 ^ 2;
	t28 = 0.1e1 / t30;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, -0.2e1 * t48 + 0.2e1 * (t28 * t46 + (t28 * t47 - t34 * t48) * t36) * t36, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:35
	% EndTime: 2019-10-09 21:44:36
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (153->13), mult. (451->28), div. (25->4), fcn. (552->7), ass. (0->20)
	t63 = sin(pkin(11));
	t64 = cos(pkin(11));
	t65 = cos(pkin(10));
	t66 = sin(qJ(2));
	t67 = cos(qJ(2));
	t78 = cos(pkin(6));
	t74 = t67 * t78;
	t75 = t66 * t78;
	t77 = sin(pkin(10));
	t54 = (t63 * t75 - t64 * t74) * t77 + t65 * (-t67 * t63 - t66 * t64);
	t48 = t54 * qJD(2);
	t72 = (-t63 * t74 - t64 * t75) * t77 - t65 * (t66 * t63 - t67 * t64);
	t51 = 0.1e1 / t72 ^ 2;
	t84 = t51 * t54 ^ 2;
	t50 = 0.1e1 / t72;
	t83 = t50 * t84;
	t82 = t51 * t72;
	t76 = t82 * t48;
	t46 = 0.1e1 + t84;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0.2e1 * (-t50 * t72 - t84) / t46 ^ 2 * (-t48 * t83 - t76) + (-0.2e1 * t76 + (t50 - t82 - 0.2e1 * t83) * t48) / t46, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:36
	% EndTime: 2019-10-09 21:44:36
	% DurationCPUTime: 0.71s
	% Computational Cost: add. (1747->58), mult. (5333->135), div. (281->12), fcn. (6885->13), ass. (0->71)
	t185 = sin(qJ(4));
	t187 = cos(qJ(4));
	t184 = cos(pkin(6));
	t179 = sin(pkin(11));
	t182 = cos(pkin(11));
	t186 = sin(qJ(2));
	t188 = cos(qJ(2));
	t200 = t188 * t179 + t186 * t182;
	t171 = t200 * t184;
	t174 = t186 * t179 - t188 * t182;
	t180 = sin(pkin(10));
	t183 = cos(pkin(10));
	t201 = -t180 * t171 - t183 * t174;
	t181 = sin(pkin(6));
	t205 = t180 * t181;
	t197 = -t185 * t201 + t187 * t205;
	t214 = t197 * qJD(4);
	t172 = t174 * qJD(2);
	t196 = t174 * t184;
	t156 = -t180 * t200 - t183 * t196;
	t169 = t174 * t181;
	t142 = atan2(t156, t169);
	t137 = sin(t142);
	t138 = cos(t142);
	t131 = t137 * t156 + t138 * t169;
	t128 = 0.1e1 / t131;
	t148 = t185 * t205 + t187 * t201;
	t144 = 0.1e1 / t148;
	t166 = 0.1e1 / t169;
	t129 = 0.1e1 / t131 ^ 2;
	t145 = 0.1e1 / t148 ^ 2;
	t167 = 0.1e1 / t169 ^ 2;
	t153 = t156 ^ 2;
	t141 = t153 * t167 + 0.1e1;
	t139 = 0.1e1 / t141;
	t195 = qJD(2) * t171;
	t149 = t180 * t172 - t183 * t195;
	t170 = t200 * t181;
	t163 = qJD(2) * t170;
	t208 = t156 * t167;
	t122 = (t149 * t166 - t163 * t208) * t139;
	t202 = -t137 * t169 + t138 * t156;
	t119 = t202 * t122 + t137 * t149 + t138 * t163;
	t213 = t119 * t128 * t129;
	t143 = t197 ^ 2;
	t134 = t143 * t145 + 0.1e1;
	t165 = t184 * t172;
	t173 = t200 * qJD(2);
	t152 = t180 * t165 - t183 * t173;
	t135 = t148 * qJD(4) + t152 * t185;
	t209 = t145 * t197;
	t136 = t152 * t187 + t214;
	t210 = t136 * t144 * t145;
	t212 = (-t135 * t209 - t143 * t210) / t134 ^ 2;
	t158 = t180 * t196 - t183 * t200;
	t211 = t129 * t158;
	t207 = t156 * t170;
	t206 = t163 * t166 * t167;
	t199 = -t144 * t185 - t187 * t209;
	t155 = -t183 * t171 + t180 * t174;
	t198 = -t155 * t166 + t167 * t207;
	t164 = t181 * t172;
	t154 = t158 ^ 2;
	t151 = t183 * t172 + t180 * t195;
	t150 = t183 * t165 + t180 * t173;
	t132 = 0.1e1 / t134;
	t126 = t154 * t129 + 0.1e1;
	t123 = t198 * t139;
	t120 = -t202 * t123 + t137 * t155 + t138 * t170;
	t118 = 0.2e1 * t198 / t141 ^ 2 * (t149 * t208 - t153 * t206) + (0.2e1 * t206 * t207 + t150 * t166 + (-t149 * t170 - t155 * t163 + t156 * t164) * t167) * t139;
	t1 = [0, t118, 0, 0, 0, 0; 0, 0.2e1 * (-t120 * t211 - t128 * t201) / t126 ^ 2 * (t151 * t211 - t154 * t213) + (t152 * t128 + (-t119 * t201 + t120 * t151) * t129 + (-0.2e1 * t120 * t213 + ((t118 * t156 - t123 * t149 - t164 + (t123 * t169 + t155) * t122) * t138 + (-t118 * t169 + t123 * t163 + t150 + (t123 * t156 - t170) * t122) * t137) * t129) * t158) / t126, 0, 0, 0, 0; 0, 0.2e1 * t199 * t158 * t212 + (-t199 * t151 + ((qJD(4) * t144 - 0.2e1 * t197 * t210) * t187 + (-t135 * t187 + (-t136 - t214) * t185) * t145) * t158) * t132, 0, -0.2e1 * t212 - 0.2e1 * (t132 * t135 * t145 - (-t132 * t210 - t145 * t212) * t197) * t197, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:36
	% EndTime: 2019-10-09 21:44:38
	% DurationCPUTime: 1.80s
	% Computational Cost: add. (5642->114), mult. (16325->235), div. (559->12), fcn. (21250->15), ass. (0->110)
	t270 = sin(pkin(11));
	t273 = cos(pkin(11));
	t278 = sin(qJ(2));
	t281 = cos(qJ(2));
	t263 = t278 * t270 - t281 * t273;
	t275 = cos(pkin(6));
	t290 = t263 * t275;
	t257 = qJD(2) * t290;
	t295 = t281 * t270 + t278 * t273;
	t262 = t295 * qJD(2);
	t271 = sin(pkin(10));
	t274 = cos(pkin(10));
	t238 = -t274 * t257 - t271 * t262;
	t260 = t295 * t275;
	t244 = t274 * t260 - t271 * t263;
	t277 = sin(qJ(4));
	t272 = sin(pkin(6));
	t312 = t272 * t277;
	t301 = t274 * t312;
	t280 = cos(qJ(4));
	t307 = qJD(4) * t280;
	t210 = -qJD(4) * t301 + t238 * t277 + t244 * t307;
	t311 = t272 * t280;
	t232 = t244 * t277 + t274 * t311;
	t230 = t232 ^ 2;
	t259 = t295 * t272;
	t251 = t259 * t277 - t275 * t280;
	t249 = 0.1e1 / t251 ^ 2;
	t226 = t230 * t249 + 0.1e1;
	t224 = 0.1e1 / t226;
	t252 = t259 * t280 + t275 * t277;
	t258 = t263 * t272;
	t256 = qJD(2) * t258;
	t228 = qJD(4) * t252 - t256 * t277;
	t248 = 0.1e1 / t251;
	t315 = t232 * t249;
	t194 = (-t210 * t248 + t228 * t315) * t224;
	t227 = atan2(-t232, t251);
	t222 = sin(t227);
	t223 = cos(t227);
	t298 = -t222 * t251 - t223 * t232;
	t190 = t194 * t298 - t222 * t210 + t223 * t228;
	t206 = -t222 * t232 + t223 * t251;
	t203 = 0.1e1 / t206;
	t204 = 0.1e1 / t206 ^ 2;
	t329 = t190 * t203 * t204;
	t296 = -t271 * t260 - t274 * t263;
	t291 = t271 * t311 - t277 * t296;
	t328 = -0.2e1 * t291 * t329;
	t243 = -t271 * t295 - t274 * t290;
	t292 = -t243 * t248 - t258 * t315;
	t327 = t277 * t292;
	t316 = t228 * t248 * t249;
	t326 = -0.2e1 * (t210 * t315 - t230 * t316) / t226 ^ 2;
	t236 = t271 * t312 + t280 * t296;
	t246 = t271 * t290 - t274 * t295;
	t276 = sin(qJ(5));
	t279 = cos(qJ(5));
	t219 = t236 * t279 - t246 * t276;
	t215 = 0.1e1 / t219;
	t216 = 0.1e1 / t219 ^ 2;
	t297 = t271 * t257 - t274 * t262;
	t213 = qJD(4) * t291 + t280 * t297;
	t261 = t263 * qJD(2);
	t289 = t275 * t262;
	t239 = t274 * t261 + t271 * t289;
	t201 = qJD(5) * t219 + t213 * t276 + t239 * t279;
	t218 = t236 * t276 + t246 * t279;
	t214 = t218 ^ 2;
	t209 = t214 * t216 + 0.1e1;
	t320 = t216 * t218;
	t306 = qJD(5) * t218;
	t202 = t213 * t279 - t239 * t276 - t306;
	t324 = t202 * t215 * t216;
	t325 = (t201 * t320 - t214 * t324) / t209 ^ 2;
	t323 = t204 * t291;
	t212 = qJD(4) * t236 + t277 * t297;
	t322 = t212 * t204;
	t321 = t215 * t276;
	t319 = t218 * t279;
	t318 = t222 * t291;
	t317 = t223 * t291;
	t314 = t246 * t277;
	t313 = t246 * t280;
	t231 = t291 ^ 2;
	t200 = t231 * t204 + 0.1e1;
	t305 = 0.2e1 * (-t231 * t329 - t291 * t322) / t200 ^ 2;
	t304 = -0.2e1 * t325;
	t302 = t218 * t324;
	t300 = -0.2e1 * t232 * t316;
	t299 = qJD(5) * t313 - t297;
	t294 = t216 * t319 - t321;
	t234 = t244 * t280 - t301;
	t293 = -t234 * t248 + t252 * t315;
	t288 = -qJD(4) * t314 + qJD(5) * t296 + t239 * t280;
	t255 = t272 * t262;
	t237 = t271 * t261 - t274 * t289;
	t229 = -qJD(4) * t251 - t256 * t280;
	t221 = t276 * t296 + t279 * t313;
	t220 = t276 * t313 - t279 * t296;
	t211 = -qJD(4) * t232 + t238 * t280;
	t207 = 0.1e1 / t209;
	t198 = 0.1e1 / t200;
	t196 = t224 * t327;
	t195 = t293 * t224;
	t192 = (-t222 * t243 - t223 * t258) * t277 + t298 * t196;
	t191 = t195 * t298 - t222 * t234 + t223 * t252;
	t189 = t293 * t326 + (t252 * t300 - t211 * t248 + (t210 * t252 + t228 * t234 + t229 * t232) * t249) * t224;
	t187 = t326 * t327 + (t292 * t307 + (-t258 * t300 - t237 * t248 + (-t210 * t258 + t228 * t243 - t232 * t255) * t249) * t277) * t224;
	t1 = [0, t187, 0, t189, 0, 0; 0, (-t192 * t323 - t203 * t314) * t305 + ((t239 * t277 + t246 * t307) * t203 + (-t322 + t328) * t192 + (-t314 * t190 + (-t258 * t307 - t187 * t232 - t196 * t210 - t255 * t277 + (-t196 * t251 - t243 * t277) * t194) * t317 + (-t243 * t307 - t187 * t251 - t196 * t228 - t237 * t277 + (t196 * t232 + t258 * t277) * t194) * t318) * t204) * t198, 0, (-t191 * t323 - t203 * t236) * t305 + (t191 * t328 + t213 * t203 + (-t236 * t190 - t191 * t212 + (-t189 * t232 - t195 * t210 + t229 + (-t195 * t251 - t234) * t194) * t317 + (-t189 * t251 - t195 * t228 - t211 + (t195 * t232 - t252) * t194) * t318) * t204) * t198, 0, 0; 0, 0.2e1 * (-t215 * t220 + t221 * t320) * t325 + (0.2e1 * t221 * t302 + t299 * t215 * t279 + t288 * t321 + (t218 * t276 * t299 - t221 * t201 - t220 * t202 - t288 * t319) * t216) * t207, 0, -t294 * t291 * t304 + (t294 * t212 - ((-qJD(5) * t215 - 0.2e1 * t302) * t279 + (t201 * t279 + (t202 - t306) * t276) * t216) * t291) * t207, t304 + 0.2e1 * (t201 * t216 * t207 + (-t207 * t324 - t216 * t325) * t218) * t218, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:44:36
	% EndTime: 2019-10-09 21:44:40
	% DurationCPUTime: 3.59s
	% Computational Cost: add. (13772->162), mult. (39160->316), div. (784->12), fcn. (51130->15), ass. (0->128)
	t373 = sin(pkin(11));
	t377 = cos(pkin(6));
	t337 = t377 * t373;
	t375 = cos(pkin(11));
	t338 = t377 * t375;
	t378 = sin(qJ(2));
	t379 = cos(qJ(2));
	t302 = t379 * t337 + t378 * t338;
	t305 = t378 * t373 - t379 * t375;
	t374 = sin(pkin(10));
	t376 = cos(pkin(10));
	t288 = t376 * t302 - t374 * t305;
	t313 = sin(qJ(4));
	t315 = cos(qJ(4));
	t311 = sin(pkin(6));
	t350 = t311 * t376;
	t275 = -t288 * t313 - t315 * t350;
	t382 = t378 * t337 - t379 * t338;
	t299 = t382 * qJD(2);
	t325 = t379 * t373 + t378 * t375;
	t304 = t325 * qJD(2);
	t281 = -t376 * t299 - t374 * t304;
	t251 = t275 * qJD(4) + t281 * t315;
	t314 = cos(qJ(5));
	t329 = -t288 * t315 + t313 * t350;
	t312 = sin(qJ(5));
	t360 = -t374 * t325 - t376 * t382;
	t383 = t360 * t312;
	t261 = -t314 * t329 - t383;
	t381 = qJD(2) * t288;
	t233 = t261 * qJD(5) + t251 * t312 - t314 * t381;
	t345 = t360 * t314;
	t259 = -t312 * t329 + t345;
	t254 = t259 ^ 2;
	t300 = t305 * t311;
	t301 = t325 * t311;
	t332 = -t301 * t315 - t377 * t313;
	t272 = -t300 * t314 - t312 * t332;
	t270 = 0.1e1 / t272 ^ 2;
	t242 = t254 * t270 + 0.1e1;
	t240 = 0.1e1 / t242;
	t292 = -t301 * t313 + t377 * t315;
	t298 = qJD(2) * t300;
	t268 = t292 * qJD(4) - t298 * t315;
	t273 = t300 * t312 - t314 * t332;
	t297 = qJD(2) * t301;
	t247 = t273 * qJD(5) + t268 * t312 - t297 * t314;
	t269 = 0.1e1 / t272;
	t364 = t259 * t270;
	t220 = (-t233 * t269 + t247 * t364) * t240;
	t243 = atan2(-t259, t272);
	t238 = sin(t243);
	t239 = cos(t243);
	t336 = -t238 * t272 - t239 * t259;
	t216 = t336 * t220 - t238 * t233 + t239 * t247;
	t232 = -t238 * t259 + t239 * t272;
	t229 = 0.1e1 / t232;
	t230 = 0.1e1 / t232 ^ 2;
	t386 = t216 * t229 * t230;
	t290 = -t376 * t325 + t374 * t382;
	t326 = -t374 * t302 - t376 * t305;
	t349 = t311 * t374;
	t328 = -t313 * t349 - t315 * t326;
	t262 = t290 * t314 - t312 * t328;
	t385 = -0.2e1 * t262;
	t344 = 0.2e1 * t262 * t386;
	t333 = -t269 * t275 + t292 * t364;
	t384 = t312 * t333;
	t367 = t247 * t269 * t270;
	t380 = -0.2e1 * (t233 * t364 - t254 * t367) / t242 ^ 2;
	t263 = -t290 * t312 - t314 * t328;
	t256 = 0.1e1 / t263;
	t257 = 0.1e1 / t263 ^ 2;
	t277 = -t313 * t326 + t315 * t349;
	t274 = t277 ^ 2;
	t366 = t257 * t274;
	t246 = 0.1e1 + t366;
	t327 = t374 * t299 - t376 * t304;
	t252 = t328 * qJD(4) - t313 * t327;
	t253 = t277 * qJD(4) + t315 * t327;
	t282 = t326 * qJD(2);
	t236 = -t262 * qJD(5) + t253 * t314 + t282 * t312;
	t370 = t236 * t256 * t257;
	t352 = t274 * t370;
	t365 = t257 * t277;
	t372 = (t252 * t365 - t352) / t246 ^ 2;
	t371 = t230 * t262;
	t369 = t238 * t262;
	t368 = t239 * t262;
	t363 = t277 * t312;
	t362 = t312 * t315;
	t361 = t314 * t315;
	t359 = qJD(4) * t313;
	t358 = qJD(5) * t312;
	t357 = qJD(5) * t314;
	t356 = qJD(5) * t315;
	t255 = t262 ^ 2;
	t228 = t230 * t255 + 0.1e1;
	t235 = t263 * qJD(5) + t253 * t312 - t282 * t314;
	t355 = 0.2e1 * (t235 * t371 - t255 * t386) / t228 ^ 2;
	t354 = 0.2e1 * t372;
	t351 = t277 * t370;
	t346 = t315 * t360;
	t343 = -0.2e1 * t259 * t367;
	t335 = -t261 * t269 + t273 * t364;
	t264 = -t288 * t314 + t312 * t346;
	t280 = -t300 * t362 - t301 * t314;
	t334 = -t264 * t269 + t280 * t364;
	t267 = t332 * qJD(4) + t298 * t313;
	t266 = t290 * t361 + t312 * t326;
	t265 = t290 * t362 - t314 * t326;
	t250 = t329 * qJD(4) - t281 * t313;
	t249 = (-t300 * t356 + t298) * t314 + (qJD(5) * t301 - t297 * t315 + t300 * t359) * t312;
	t248 = -t272 * qJD(5) + t268 * t314 + t297 * t312;
	t244 = 0.1e1 / t246;
	t237 = -t281 * t314 + t288 * t358 + t346 * t357 - t359 * t383 - t362 * t381;
	t234 = -qJD(5) * t345 + t251 * t314 + t312 * t381 + t329 * t358;
	t226 = 0.1e1 / t228;
	t225 = t240 * t384;
	t224 = t334 * t240;
	t222 = t335 * t240;
	t219 = (-t238 * t275 + t239 * t292) * t312 + t336 * t225;
	t218 = t336 * t224 - t238 * t264 + t239 * t280;
	t217 = t336 * t222 - t238 * t261 + t239 * t273;
	t215 = t334 * t380 + (t280 * t343 - t237 * t269 + (t233 * t280 + t247 * t264 + t249 * t259) * t270) * t240;
	t213 = t335 * t380 + (t273 * t343 - t234 * t269 + (t233 * t273 + t247 * t261 + t248 * t259) * t270) * t240;
	t212 = t380 * t384 + (t333 * t357 + (t292 * t343 - t250 * t269 + (t233 * t292 + t247 * t275 + t259 * t267) * t270) * t312) * t240;
	t1 = [0, t215, 0, t212, t213, 0; 0, (t218 * t371 - t229 * t265) * t355 + (t218 * t344 + (-t265 * t216 - t218 * t235 - (-t215 * t259 - t224 * t233 + t249 + (-t224 * t272 - t264) * t220) * t368 - (-t215 * t272 - t224 * t247 - t237 + (t224 * t259 - t280) * t220) * t369) * t230 + ((t290 * t356 - t327) * t314 + (qJD(5) * t326 - t282 * t315 - t290 * t359) * t312) * t229) * t226, 0, (t219 * t371 - t229 * t363) * t355 + ((t252 * t312 + t277 * t357) * t229 + t219 * t344 + (-t219 * t235 - t363 * t216 - (t292 * t357 - t212 * t259 - t225 * t233 + t267 * t312 + (-t225 * t272 - t275 * t312) * t220) * t368 - (-t275 * t357 - t212 * t272 - t225 * t247 - t250 * t312 + (t225 * t259 - t292 * t312) * t220) * t369) * t230) * t226, (t217 * t371 - t229 * t263) * t355 + (t217 * t344 + t236 * t229 + (-t263 * t216 - t217 * t235 - (-t213 * t259 - t222 * t233 + t248 + (-t222 * t272 - t261) * t220) * t368 - (-t213 * t272 - t222 * t247 - t234 + (t222 * t259 - t273) * t220) * t369) * t230) * t226, 0; 0, (t256 * t290 * t313 + t266 * t365) * t354 + (0.2e1 * t266 * t351 + (-qJD(4) * t290 * t315 + t282 * t313) * t256 + (-(-t282 * t361 + t312 * t327 + t326 * t357) * t277 - t266 * t252 + (t313 * t236 - (-t312 * t356 - t314 * t359) * t277) * t290) * t257) * t244, 0, (-t256 * t328 + t314 * t366) * t354 + (0.2e1 * t314 * t352 - t253 * t256 + (-0.2e1 * t252 * t277 * t314 - t236 * t328 + t274 * t358) * t257) * t244, t365 * t372 * t385 + (t351 * t385 + (t235 * t277 + t252 * t262) * t257) * t244, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end