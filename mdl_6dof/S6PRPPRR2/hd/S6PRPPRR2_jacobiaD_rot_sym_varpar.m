% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPPRR2
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
%   Wie in S6PRPPRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:26
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPPRR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
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
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
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
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (1331->44), mult. (4118->102), div. (242->14), fcn. (5370->11), ass. (0->57)
	t122 = cos(pkin(6));
	t117 = sin(pkin(11));
	t120 = cos(pkin(11));
	t123 = sin(qJ(2));
	t124 = cos(qJ(2));
	t136 = t124 * t117 + t123 * t120;
	t110 = t136 * t122;
	t113 = t123 * t117 - t124 * t120;
	t111 = t113 * qJD(2);
	t119 = sin(pkin(6));
	t108 = t113 * t119;
	t118 = sin(pkin(10));
	t121 = cos(pkin(10));
	t134 = t113 * t122;
	t95 = -t118 * t136 - t121 * t134;
	t87 = atan2(t95, t108);
	t82 = sin(t87);
	t83 = cos(t87);
	t81 = t83 * t108 + t82 * t95;
	t78 = 0.1e1 / t81;
	t97 = t118 * t134 - t121 * t136;
	t149 = -0.2e1 * t97;
	t105 = 0.1e1 / t108;
	t106 = 0.1e1 / t108 ^ 2;
	t79 = 0.1e1 / t81 ^ 2;
	t109 = t136 * t119;
	t102 = qJD(2) * t109;
	t139 = -t108 * t82 + t83 * t95;
	t146 = t106 * t95;
	t92 = t95 ^ 2;
	t86 = t92 * t106 + 0.1e1;
	t84 = 0.1e1 / t86;
	t133 = qJD(2) * t110;
	t89 = t118 * t111 - t121 * t133;
	t72 = (-t102 * t146 + t105 * t89) * t84;
	t70 = t83 * t102 + t139 * t72 + t82 * t89;
	t148 = t70 * t78 * t79;
	t147 = t79 * t97;
	t145 = t109 * t95;
	t144 = t102 * t105 * t106;
	t140 = 0.1e1 / t118 ^ 2 / t119 ^ 2;
	t104 = t122 * t111;
	t112 = t136 * qJD(2);
	t138 = t118 * t104 - t121 * t112;
	t137 = -t118 * t110 - t121 * t113;
	t94 = -t121 * t110 + t118 * t113;
	t135 = -t105 * t94 + t106 * t145;
	t103 = t119 * t111;
	t93 = t97 ^ 2;
	t91 = t121 * t111 + t118 * t133;
	t90 = t121 * t104 + t118 * t112;
	t88 = t137 ^ 2 * t140 + 0.1e1;
	t76 = t93 * t79 + 0.1e1;
	t73 = t135 * t84;
	t71 = t83 * t109 - t139 * t73 + t82 * t94;
	t69 = 0.2e1 * t135 / t86 ^ 2 * (-t92 * t144 + t89 * t146) + (0.2e1 * t144 * t145 + t105 * t90 + (-t102 * t94 + t103 * t95 - t109 * t89) * t106) * t84;
	t1 = [0, t69, 0, 0, 0, 0; 0, 0.2e1 * (-t137 * t78 - t71 * t147) / t76 ^ 2 * (t91 * t147 - t93 * t148) + (t138 * t78 + t71 * t148 * t149 + (-t137 * t70 + t71 * t91 + ((t69 * t95 - t73 * t89 - t103 + (t108 * t73 + t94) * t72) * t83 + (t102 * t73 - t108 * t69 + t90 + (t73 * t95 - t109) * t72) * t82) * t97) * t79) / t76, 0, 0, 0, 0; 0, (t91 / t88 + 0.1e1 / t88 ^ 2 * t137 * t138 * t140 * t149) / t118 / t119, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:20
	% DurationCPUTime: 0.70s
	% Computational Cost: add. (1747->59), mult. (5333->135), div. (281->12), fcn. (6885->13), ass. (0->71)
	t181 = sin(pkin(11));
	t184 = cos(pkin(11));
	t188 = sin(qJ(2));
	t190 = cos(qJ(2));
	t177 = t188 * t181 - t190 * t184;
	t186 = cos(pkin(6));
	t198 = t177 * t186;
	t201 = t190 * t181 + t188 * t184;
	t174 = t201 * t186;
	t182 = sin(pkin(10));
	t185 = cos(pkin(10));
	t159 = t185 * t174 - t182 * t177;
	t183 = sin(pkin(6));
	t173 = t201 * t183;
	t143 = atan2(-t159, t173);
	t138 = sin(t143);
	t139 = cos(t143);
	t132 = -t138 * t159 + t139 * t173;
	t129 = 0.1e1 / t132;
	t162 = t182 * t198 - t185 * t201;
	t187 = sin(qJ(5));
	t189 = cos(qJ(5));
	t208 = t182 * t183;
	t149 = -t162 * t187 + t189 * t208;
	t145 = 0.1e1 / t149;
	t169 = 0.1e1 / t173;
	t130 = 0.1e1 / t132 ^ 2;
	t146 = 0.1e1 / t149 ^ 2;
	t170 = 0.1e1 / t173 ^ 2;
	t156 = t159 ^ 2;
	t142 = t156 * t170 + 0.1e1;
	t140 = 0.1e1 / t142;
	t168 = qJD(2) * t198;
	t176 = t201 * qJD(2);
	t151 = -t185 * t168 - t182 * t176;
	t172 = t177 * t183;
	t167 = qJD(2) * t172;
	t211 = t159 * t170;
	t123 = (-t151 * t169 - t167 * t211) * t140;
	t203 = -t138 * t173 - t139 * t159;
	t120 = t203 * t123 - t138 * t151 - t139 * t167;
	t216 = t120 * t129 * t130;
	t148 = t162 * t189 + t187 * t208;
	t144 = t148 ^ 2;
	t135 = t144 * t146 + 0.1e1;
	t175 = t177 * qJD(2);
	t197 = t186 * t176;
	t152 = t185 * t175 + t182 * t197;
	t137 = t149 * qJD(5) + t152 * t189;
	t212 = t146 * t148;
	t204 = qJD(5) * t148;
	t136 = -t152 * t187 - t204;
	t213 = t136 * t145 * t146;
	t215 = (t137 * t212 - t144 * t213) / t135 ^ 2;
	t202 = -t182 * t174 - t185 * t177;
	t214 = t130 * t202;
	t210 = t159 * t172;
	t209 = t167 * t169 * t170;
	t153 = t182 * t168 - t185 * t176;
	t200 = t145 * t189 + t187 * t212;
	t158 = -t182 * t201 - t185 * t198;
	t199 = -t158 * t169 - t170 * t210;
	t166 = t183 * t176;
	t157 = t202 ^ 2;
	t150 = t182 * t175 - t185 * t197;
	t133 = 0.1e1 / t135;
	t127 = t157 * t130 + 0.1e1;
	t124 = t199 * t140;
	t121 = t203 * t124 - t138 * t158 - t139 * t172;
	t119 = -0.2e1 * t199 / t142 ^ 2 * (t151 * t211 + t156 * t209) + (-0.2e1 * t209 * t210 - t150 * t169 + (-t151 * t172 - t158 * t167 - t159 * t166) * t170) * t140;
	t1 = [0, t119, 0, 0, 0, 0; 0, 0.2e1 * (t121 * t214 - t129 * t162) / t127 ^ 2 * (t153 * t214 - t157 * t216) + (t152 * t129 + (-t162 * t120 - t121 * t153) * t130 + (0.2e1 * t121 * t216 + (-(-t119 * t159 - t124 * t151 - t166 + (-t124 * t173 - t158) * t123) * t139 - (-t119 * t173 + t124 * t167 - t150 + (t124 * t159 + t172) * t123) * t138) * t130) * t202) / t127, 0, 0, 0, 0; 0, 0.2e1 * t200 * t202 * t215 + (-t200 * t153 + ((qJD(5) * t145 + 0.2e1 * t148 * t213) * t187 + (-t137 * t187 + (t136 - t204) * t189) * t146) * t202) * t133, 0, 0, -0.2e1 * t215 + 0.2e1 * (t133 * t137 * t146 + (-t133 * t213 - t146 * t215) * t148) * t148, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:21
	% DurationCPUTime: 1.77s
	% Computational Cost: add. (5642->114), mult. (16325->239), div. (559->12), fcn. (21250->15), ass. (0->110)
	t277 = cos(pkin(6));
	t272 = sin(pkin(11));
	t275 = cos(pkin(11));
	t280 = sin(qJ(2));
	t283 = cos(qJ(2));
	t296 = t283 * t272 + t280 * t275;
	t265 = t296 * t277;
	t261 = qJD(2) * t265;
	t269 = t280 * t272 - t283 * t275;
	t267 = t269 * qJD(2);
	t273 = sin(pkin(10));
	t276 = cos(pkin(10));
	t240 = t276 * t261 - t273 * t267;
	t274 = sin(pkin(6));
	t282 = cos(qJ(5));
	t279 = sin(qJ(5));
	t264 = t269 * t277;
	t291 = t276 * t264 + t273 * t296;
	t290 = t291 * t279;
	t216 = -qJD(5) * t290 + (t276 * t274 * qJD(5) + t240) * t282;
	t312 = t274 * t279;
	t237 = t276 * t312 + t291 * t282;
	t234 = t237 ^ 2;
	t262 = t269 * t274;
	t298 = t262 * t282 - t277 * t279;
	t254 = 0.1e1 / t298 ^ 2;
	t229 = t234 * t254 + 0.1e1;
	t227 = 0.1e1 / t229;
	t257 = t262 * t279 + t277 * t282;
	t263 = t296 * t274;
	t258 = qJD(2) * t263;
	t232 = t257 * qJD(5) - t258 * t282;
	t253 = 0.1e1 / t298;
	t315 = t237 * t254;
	t197 = (-t216 * t253 - t232 * t315) * t227;
	t230 = atan2(t237, -t298);
	t225 = sin(t230);
	t226 = cos(t230);
	t300 = t225 * t298 + t226 * t237;
	t193 = t300 * t197 + t225 * t216 + t226 * t232;
	t209 = t225 * t237 - t226 * t298;
	t206 = 0.1e1 / t209;
	t207 = 0.1e1 / t209 ^ 2;
	t329 = t193 * t206 * t207;
	t251 = t273 * t264 - t276 * t296;
	t235 = t251 * t282 + t273 * t312;
	t328 = 0.2e1 * t235 * t329;
	t316 = t232 * t253 * t254;
	t327 = (t216 * t315 + t234 * t316) / t229 ^ 2;
	t249 = t276 * t265 - t273 * t269;
	t293 = -t249 * t253 + t263 * t315;
	t326 = t282 * t293;
	t311 = t274 * t282;
	t236 = -t251 * t279 + t273 * t311;
	t278 = sin(qJ(6));
	t281 = cos(qJ(6));
	t297 = -t273 * t265 - t276 * t269;
	t222 = t236 * t281 + t278 * t297;
	t218 = 0.1e1 / t222;
	t219 = 0.1e1 / t222 ^ 2;
	t242 = t273 * t261 + t276 * t267;
	t213 = -t235 * qJD(5) - t242 * t279;
	t260 = t277 * t267;
	t268 = t296 * qJD(2);
	t299 = t273 * t260 - t276 * t268;
	t204 = t222 * qJD(6) + t213 * t278 - t281 * t299;
	t221 = t236 * t278 - t281 * t297;
	t217 = t221 ^ 2;
	t212 = t217 * t219 + 0.1e1;
	t320 = t219 * t221;
	t307 = qJD(6) * t221;
	t205 = t213 * t281 + t278 * t299 - t307;
	t324 = t205 * t218 * t219;
	t325 = (t204 * t320 - t217 * t324) / t212 ^ 2;
	t323 = t207 * t235;
	t214 = t236 * qJD(5) + t242 * t282;
	t322 = t214 * t207;
	t321 = t218 * t278;
	t319 = t221 * t281;
	t318 = t225 * t235;
	t317 = t226 * t235;
	t314 = t297 * t279;
	t313 = t297 * t282;
	t308 = qJD(5) * t279;
	t233 = t235 ^ 2;
	t203 = t233 * t207 + 0.1e1;
	t306 = 0.2e1 * (-t233 * t329 + t235 * t322) / t203 ^ 2;
	t305 = -0.2e1 * t325;
	t303 = t221 * t324;
	t302 = t237 * t316;
	t301 = qJD(6) * t314 - t242;
	t295 = t219 * t319 - t321;
	t238 = t276 * t311 - t290;
	t294 = t238 * t253 + t257 * t315;
	t292 = qJD(5) * t313 + qJD(6) * t251 + t279 * t299;
	t259 = t274 * t267;
	t241 = -t276 * t260 - t273 * t268;
	t231 = t298 * qJD(5) + t258 * t279;
	t224 = t251 * t278 + t281 * t314;
	t223 = -t251 * t281 + t278 * t314;
	t215 = t237 * qJD(5) + t240 * t279;
	t210 = 0.1e1 / t212;
	t201 = 0.1e1 / t203;
	t199 = t227 * t326;
	t198 = t294 * t227;
	t195 = (t225 * t249 - t226 * t263) * t282 + t300 * t199;
	t194 = -t300 * t198 + t225 * t238 + t226 * t257;
	t192 = 0.2e1 * t294 * t327 + (-0.2e1 * t257 * t302 + t215 * t253 + (-t216 * t257 - t231 * t237 - t232 * t238) * t254) * t227;
	t190 = -0.2e1 * t326 * t327 + (-t293 * t308 + (0.2e1 * t263 * t302 - t241 * t253 + (t216 * t263 - t232 * t249 - t237 * t259) * t254) * t282) * t227;
	t1 = [0, t190, 0, 0, t192, 0; 0, (t195 * t323 + t206 * t313) * t306 + ((-t282 * t299 + t297 * t308) * t206 + (-t322 + t328) * t195 + (t313 * t193 - (t263 * t308 + t190 * t237 + t199 * t216 + t259 * t282 + (t199 * t298 + t249 * t282) * t197) * t317 - (-t249 * t308 + t190 * t298 - t199 * t232 + t241 * t282 + (-t199 * t237 + t263 * t282) * t197) * t318) * t207) * t201, 0, 0, (t194 * t323 - t206 * t236) * t306 + (t194 * t328 + t213 * t206 + (-t236 * t193 - t194 * t214 - (t192 * t237 - t198 * t216 + t231 + (-t198 * t298 + t238) * t197) * t317 - (t192 * t298 + t198 * t232 - t215 + (t198 * t237 - t257) * t197) * t318) * t207) * t201, 0; 0, 0.2e1 * (-t218 * t223 + t224 * t320) * t325 + (0.2e1 * t224 * t303 + t301 * t218 * t281 + t292 * t321 + (t301 * t221 * t278 - t224 * t204 - t223 * t205 - t292 * t319) * t219) * t210, 0, 0, t295 * t235 * t305 + (t295 * t214 + ((-qJD(6) * t218 - 0.2e1 * t303) * t281 + (t204 * t281 + (t205 - t307) * t278) * t219) * t235) * t210, t305 + 0.2e1 * (t204 * t219 * t210 + (-t210 * t324 - t219 * t325) * t221) * t221;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end