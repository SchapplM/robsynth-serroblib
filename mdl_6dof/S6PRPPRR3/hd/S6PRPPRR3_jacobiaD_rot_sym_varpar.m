% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPPRR3
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
%   Wie in S6PRPPRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:28
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPPRR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:09
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (543->41), mult. (1611->102), div. (383->14), fcn. (2156->9), ass. (0->52)
	t119 = cos(pkin(6));
	t94 = sin(pkin(10));
	t113 = t94 * t119;
	t118 = cos(pkin(10));
	t96 = sin(qJ(2));
	t97 = cos(qJ(2));
	t87 = -t96 * t113 + t118 * t97;
	t82 = 0.1e1 / t87 ^ 2;
	t95 = sin(pkin(6));
	t105 = t95 ^ 2;
	t126 = t82 * t94 ^ 2 * t105;
	t120 = t95 * t97;
	t112 = t119 * t118;
	t83 = -t97 * t112 + t94 * t96;
	t70 = atan2(-t83, -t120);
	t68 = sin(t70);
	t69 = cos(t70);
	t66 = -t69 * t120 - t68 * t83;
	t63 = 0.1e1 / t66;
	t91 = 0.1e1 / t97;
	t64 = 0.1e1 / t66 ^ 2;
	t92 = 0.1e1 / t97 ^ 2;
	t111 = t68 * t120 - t69 * t83;
	t116 = t69 * t95 * t96;
	t122 = t92 * t96;
	t115 = t83 * t122;
	t80 = t83 ^ 2;
	t90 = 0.1e1 / t105;
	t75 = t80 * t90 * t92 + 0.1e1;
	t72 = 0.1e1 / t75;
	t89 = 0.1e1 / t95;
	t123 = t72 * t89;
	t85 = t96 * t112 + t94 * t97;
	t77 = t85 * qJD(2);
	t58 = (qJD(2) * t115 + t77 * t91) * t123;
	t56 = qJD(2) * t116 + t111 * t58 - t68 * t77;
	t125 = t56 * t63 * t64;
	t109 = -t97 * t113 - t118 * t96;
	t124 = t64 * t109;
	t78 = t109 * qJD(2);
	t114 = t109 / t87 * t82 * t78;
	t110 = t85 * t91 + t115;
	t93 = t91 * t92;
	t81 = t109 ^ 2;
	t79 = t87 * qJD(2);
	t76 = qJD(2) * t83;
	t74 = 0.1e1 + t126;
	t62 = t81 * t64 + 0.1e1;
	t59 = t110 * t123;
	t57 = t111 * t59 - t68 * t85 + t116;
	t55 = (-0.2e1 * t110 / t75 ^ 2 * (qJD(2) * t80 * t93 * t96 + t77 * t83 * t92) * t90 + (t77 * t122 - t76 * t91 + (t85 * t122 + (0.2e1 * t93 * t96 ^ 2 + t91) * t83) * qJD(2)) * t72) * t89;
	t1 = [0, t55, 0, 0, 0, 0; 0, 0.2e1 * (-t57 * t124 - t63 * t87) / t62 ^ 2 * (-t79 * t124 - t81 * t125) + (t78 * t63 + (-t87 * t56 - t57 * t79) * t64 - (0.2e1 * t57 * t125 + (-(qJD(2) * t120 - t55 * t83 - t59 * t77 + (t59 * t120 - t85) * t58) * t69 - (t58 * t59 * t83 + t76 + (t55 * t97 + (-qJD(2) * t59 - t58) * t96) * t95) * t68) * t64) * t109) / t62, 0, 0, 0, 0; 0, (0.2e1 / t74 ^ 2 * t114 * t126 + (-t79 * t82 - 0.2e1 * t114) / t74) * t94 * t95, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (153->13), mult. (451->29), div. (25->4), fcn. (552->7), ass. (0->21)
	t73 = cos(pkin(10));
	t75 = sin(qJ(2));
	t76 = cos(qJ(2));
	t82 = sin(pkin(10)) * cos(pkin(6));
	t67 = t73 * t75 + t76 * t82;
	t68 = t73 * t76 - t75 * t82;
	t70 = sin(pkin(11));
	t72 = cos(pkin(11));
	t62 = t67 * t70 + t68 * t72;
	t57 = 0.1e1 / t62 ^ 2;
	t79 = -t67 * t72 + t68 * t70;
	t89 = t57 * t79 ^ 2;
	t65 = t67 * qJD(2);
	t66 = t68 * qJD(2);
	t54 = -t65 * t72 + t66 * t70;
	t56 = 0.1e1 / t62;
	t86 = t54 * t56;
	t88 = t86 * t89;
	t85 = t79 * (-t65 * t70 - t66 * t72);
	t52 = 0.1e1 + t89;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0.2e1 * (t56 * t62 + t89) / t52 ^ 2 * (t57 * t85 - t88) + (-t86 + 0.2e1 * t88 + (t54 * t62 - 0.2e1 * t85) * t57) / t52, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:09
	% DurationCPUTime: 0.79s
	% Computational Cost: add. (1747->63), mult. (5333->142), div. (281->12), fcn. (6885->13), ass. (0->74)
	t183 = sin(pkin(11));
	t184 = sin(pkin(10));
	t188 = sin(qJ(2));
	t190 = cos(qJ(2));
	t220 = cos(pkin(10));
	t221 = cos(pkin(6));
	t205 = t221 * t220;
	t198 = -t184 * t190 - t188 * t205;
	t222 = t198 * t183;
	t185 = sin(pkin(6));
	t186 = cos(pkin(11));
	t175 = (t183 * t188 + t186 * t190) * t185;
	t209 = t184 * t221;
	t181 = -t188 * t209 + t220 * t190;
	t200 = t184 * t188 - t190 * t205;
	t159 = -t200 * t186 - t222;
	t144 = atan2(-t159, t175);
	t139 = sin(t144);
	t140 = cos(t144);
	t133 = -t139 * t159 + t140 * t175;
	t130 = 0.1e1 / t133;
	t199 = -t220 * t188 - t190 * t209;
	t164 = t181 * t186 - t183 * t199;
	t187 = sin(qJ(5));
	t189 = cos(qJ(5));
	t211 = t184 * t185;
	t150 = t164 * t189 - t187 * t211;
	t146 = 0.1e1 / t150;
	t172 = 0.1e1 / t175;
	t131 = 0.1e1 / t133 ^ 2;
	t147 = 0.1e1 / t150 ^ 2;
	t173 = 0.1e1 / t175 ^ 2;
	t156 = t159 ^ 2;
	t143 = t156 * t173 + 0.1e1;
	t141 = 0.1e1 / t143;
	t177 = t200 * qJD(2);
	t197 = t198 * t186;
	t151 = qJD(2) * t197 - t177 * t183;
	t176 = (t183 * t190 - t186 * t188) * t185;
	t168 = qJD(2) * t176;
	t214 = t159 * t173;
	t124 = (-t151 * t172 + t168 * t214) * t141;
	t204 = -t139 * t175 - t140 * t159;
	t121 = t204 * t124 - t139 * t151 + t140 * t168;
	t219 = t121 * t130 * t131;
	t149 = t164 * t187 + t189 * t211;
	t145 = t149 ^ 2;
	t136 = t145 * t147 + 0.1e1;
	t178 = t199 * qJD(2);
	t179 = t181 * qJD(2);
	t154 = t178 * t186 + t179 * t183;
	t137 = t150 * qJD(5) + t154 * t187;
	t215 = t147 * t149;
	t210 = qJD(5) * t149;
	t138 = t154 * t189 - t210;
	t216 = t138 * t146 * t147;
	t218 = (t137 * t215 - t145 * t216) / t136 ^ 2;
	t207 = t181 * t183 + t186 * t199;
	t217 = t131 * t207;
	t213 = t159 * t176;
	t212 = t168 * t172 * t173;
	t153 = t178 * t183 - t179 * t186;
	t202 = -t146 * t187 + t189 * t215;
	t158 = -t200 * t183 + t197;
	t201 = -t158 * t172 + t173 * t213;
	t167 = qJD(2) * t175;
	t157 = t207 ^ 2;
	t152 = qJD(2) * t222 + t177 * t186;
	t134 = 0.1e1 / t136;
	t128 = t157 * t131 + 0.1e1;
	t125 = t201 * t141;
	t122 = t204 * t125 - t139 * t158 + t140 * t176;
	t120 = -0.2e1 * t201 / t143 ^ 2 * (t151 * t214 - t156 * t212) + (-0.2e1 * t212 * t213 - t152 * t172 + (t151 * t176 + t158 * t168 - t159 * t167) * t173) * t141;
	t1 = [0, t120, 0, 0, 0, 0; 0, 0.2e1 * (t122 * t217 + t130 * t164) / t128 ^ 2 * (t153 * t217 - t157 * t219) + (-t154 * t130 + (t121 * t164 - t122 * t153) * t131 + (0.2e1 * t122 * t219 + (-(-t120 * t159 - t125 * t151 - t167 + (-t125 * t175 - t158) * t124) * t140 - (-t120 * t175 - t125 * t168 - t152 + (t125 * t159 - t176) * t124) * t139) * t131) * t207) / t128, 0, 0, 0, 0; 0, 0.2e1 * t202 * t207 * t218 + (-t202 * t153 + ((qJD(5) * t146 + 0.2e1 * t149 * t216) * t189 + (-t137 * t189 + (-t138 + t210) * t187) * t147) * t207) * t134, 0, 0, -0.2e1 * t218 + 0.2e1 * (t134 * t137 * t147 + (-t134 * t216 - t147 * t218) * t149) * t149, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:09
	% EndTime: 2019-10-09 21:28:11
	% DurationCPUTime: 1.80s
	% Computational Cost: add. (5642->117), mult. (16325->241), div. (559->12), fcn. (21250->15), ass. (0->114)
	t283 = cos(pkin(6));
	t286 = sin(qJ(2));
	t336 = cos(pkin(10));
	t310 = t336 * t286;
	t280 = sin(pkin(10));
	t289 = cos(qJ(2));
	t320 = t280 * t289;
	t273 = t283 * t310 + t320;
	t279 = sin(pkin(11));
	t282 = cos(pkin(11));
	t309 = t336 * t289;
	t321 = t280 * t286;
	t297 = t283 * t309 - t321;
	t252 = t273 * t282 - t279 * t297;
	t285 = sin(qJ(5));
	t288 = cos(qJ(5));
	t281 = sin(pkin(6));
	t311 = t281 * t336;
	t241 = t252 * t288 + t285 * t311;
	t268 = t297 * qJD(2);
	t269 = t273 * qJD(2);
	t244 = t268 * t282 + t269 * t279;
	t217 = qJD(5) * t241 + t244 * t285;
	t298 = -t252 * t285 + t288 * t311;
	t237 = t298 ^ 2;
	t267 = (t279 * t289 - t282 * t286) * t281;
	t260 = -t267 * t285 + t283 * t288;
	t258 = 0.1e1 / t260 ^ 2;
	t233 = t237 * t258 + 0.1e1;
	t231 = 0.1e1 / t233;
	t261 = -t267 * t288 - t283 * t285;
	t266 = (t279 * t286 + t282 * t289) * t281;
	t264 = qJD(2) * t266;
	t235 = qJD(5) * t261 + t264 * t285;
	t257 = 0.1e1 / t260;
	t325 = t298 * t258;
	t201 = (-t217 * t257 - t235 * t325) * t231;
	t234 = atan2(t298, t260);
	t229 = sin(t234);
	t230 = cos(t234);
	t304 = -t229 * t260 + t230 * t298;
	t197 = t201 * t304 - t229 * t217 + t230 * t235;
	t213 = t229 * t298 + t230 * t260;
	t210 = 0.1e1 / t213;
	t211 = 0.1e1 / t213 ^ 2;
	t340 = t197 * t210 * t211;
	t275 = -t283 * t321 + t309;
	t299 = -t283 * t320 - t310;
	t256 = t275 * t282 - t279 * t299;
	t322 = t280 * t281;
	t242 = t256 * t285 + t288 * t322;
	t339 = 0.2e1 * t242 * t340;
	t251 = t273 * t279 + t282 * t297;
	t300 = -t251 * t257 - t266 * t325;
	t338 = t285 * t300;
	t326 = t235 * t257 * t258;
	t337 = -0.2e1 * (-t217 * t325 - t237 * t326) / t233 ^ 2;
	t313 = t285 * t322;
	t243 = t256 * t288 - t313;
	t284 = sin(qJ(6));
	t287 = cos(qJ(6));
	t306 = t275 * t279 + t282 * t299;
	t226 = t243 * t287 + t284 * t306;
	t222 = 0.1e1 / t226;
	t223 = 0.1e1 / t226 ^ 2;
	t270 = t299 * qJD(2);
	t271 = t275 * qJD(2);
	t247 = t270 * t282 + t271 * t279;
	t220 = -qJD(5) * t242 + t247 * t288;
	t307 = t270 * t279 - t271 * t282;
	t208 = qJD(6) * t226 + t220 * t284 - t287 * t307;
	t225 = t243 * t284 - t287 * t306;
	t221 = t225 ^ 2;
	t216 = t221 * t223 + 0.1e1;
	t330 = t223 * t225;
	t318 = qJD(6) * t225;
	t209 = t220 * t287 + t284 * t307 - t318;
	t334 = t209 * t222 * t223;
	t335 = (t208 * t330 - t221 * t334) / t216 ^ 2;
	t333 = t211 * t242;
	t319 = qJD(5) * t288;
	t219 = -qJD(5) * t313 + t247 * t285 + t256 * t319;
	t332 = t219 * t211;
	t331 = t222 * t284;
	t329 = t225 * t287;
	t328 = t229 * t242;
	t327 = t230 * t242;
	t324 = t306 * t285;
	t323 = t306 * t288;
	t238 = t242 ^ 2;
	t207 = t238 * t211 + 0.1e1;
	t317 = 0.2e1 * (-t238 * t340 + t242 * t332) / t207 ^ 2;
	t316 = -0.2e1 * t335;
	t314 = t225 * t334;
	t308 = 0.2e1 * t298 * t326;
	t305 = qJD(6) * t323 + t247;
	t302 = t329 * t223 - t331;
	t301 = -t241 * t257 - t261 * t325;
	t296 = -qJD(5) * t324 - qJD(6) * t256 + t288 * t307;
	t263 = qJD(2) * t267;
	t245 = t268 * t279 - t269 * t282;
	t236 = -qJD(5) * t260 + t264 * t288;
	t228 = -t256 * t284 + t287 * t323;
	t227 = t256 * t287 + t284 * t323;
	t218 = qJD(5) * t298 + t244 * t288;
	t214 = 0.1e1 / t216;
	t205 = 0.1e1 / t207;
	t203 = t231 * t338;
	t202 = t301 * t231;
	t199 = (-t229 * t251 + t230 * t266) * t285 + t304 * t203;
	t198 = t202 * t304 - t229 * t241 + t230 * t261;
	t196 = t301 * t337 + (t261 * t308 - t218 * t257 + (t217 * t261 + t235 * t241 - t236 * t298) * t258) * t231;
	t194 = t337 * t338 + (t300 * t319 + (t266 * t308 - t245 * t257 + (t217 * t266 + t235 * t251 - t263 * t298) * t258) * t285) * t231;
	t1 = [0, t194, 0, 0, t196, 0; 0, (t199 * t333 - t210 * t324) * t317 + ((t285 * t307 + t306 * t319) * t210 + (-t332 + t339) * t199 + (-t324 * t197 - (t266 * t319 + t194 * t298 - t203 * t217 + t263 * t285 + (-t203 * t260 - t251 * t285) * t201) * t327 - (-t251 * t319 - t194 * t260 - t203 * t235 - t245 * t285 + (-t203 * t298 - t266 * t285) * t201) * t328) * t211) * t205, 0, 0, (t198 * t333 - t210 * t243) * t317 + (t198 * t339 + t220 * t210 + (-t243 * t197 - t198 * t219 - (t196 * t298 - t202 * t217 + t236 + (-t202 * t260 - t241) * t201) * t327 - (-t196 * t260 - t202 * t235 - t218 + (-t202 * t298 - t261) * t201) * t328) * t211) * t205, 0; 0, 0.2e1 * (-t222 * t227 + t228 * t330) * t335 + (0.2e1 * t228 * t314 + t305 * t222 * t287 + t296 * t331 + (t225 * t284 * t305 - t228 * t208 - t227 * t209 - t296 * t329) * t223) * t214, 0, 0, t302 * t242 * t316 + (t302 * t219 + ((-qJD(6) * t222 - 0.2e1 * t314) * t287 + (t208 * t287 + (t209 - t318) * t284) * t223) * t242) * t214, t316 + 0.2e1 * (t208 * t223 * t214 + (-t214 * t334 - t223 * t335) * t225) * t225;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end