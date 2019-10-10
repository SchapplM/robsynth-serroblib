% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR5
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
%   Wie in S6PRPRRR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:01
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(11));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(11)) * cos(pkin(6));
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
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:17
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (536->40), mult. (1574->99), div. (384->14), fcn. (2127->9), ass. (0->50)
	t90 = sin(pkin(6));
	t85 = 0.1e1 / t90 ^ 2;
	t89 = sin(pkin(11));
	t118 = 0.1e1 / t89 ^ 2 * t85;
	t92 = cos(qJ(2));
	t110 = t90 * t92;
	t108 = cos(pkin(11));
	t109 = cos(pkin(6));
	t104 = t109 * t108;
	t91 = sin(qJ(2));
	t78 = -t92 * t104 + t89 * t91;
	t67 = atan2(-t78, -t110);
	t65 = sin(t67);
	t66 = cos(t67);
	t63 = -t66 * t110 - t65 * t78;
	t60 = 0.1e1 / t63;
	t84 = 0.1e1 / t90;
	t86 = 0.1e1 / t92;
	t61 = 0.1e1 / t63 ^ 2;
	t87 = 0.1e1 / t92 ^ 2;
	t105 = t89 * t109;
	t101 = -t92 * t105 - t108 * t91;
	t115 = -0.2e1 * t101;
	t103 = t65 * t110 - t66 * t78;
	t107 = t66 * t90 * t91;
	t111 = t87 * t91;
	t106 = t78 * t111;
	t76 = t78 ^ 2;
	t71 = t76 * t85 * t87 + 0.1e1;
	t68 = 0.1e1 / t71;
	t112 = t68 * t84;
	t80 = t91 * t104 + t89 * t92;
	t73 = t80 * qJD(2);
	t55 = (qJD(2) * t106 + t73 * t86) * t112;
	t53 = qJD(2) * t107 + t103 * t55 - t65 * t73;
	t114 = t53 * t60 * t61;
	t113 = t61 * t101;
	t102 = t80 * t86 + t106;
	t82 = -t91 * t105 + t108 * t92;
	t88 = t86 * t87;
	t77 = t101 ^ 2;
	t75 = t82 * qJD(2);
	t74 = t101 * qJD(2);
	t72 = qJD(2) * t78;
	t70 = t82 ^ 2 * t118 + 0.1e1;
	t59 = t77 * t61 + 0.1e1;
	t56 = t102 * t112;
	t54 = t103 * t56 - t65 * t80 + t107;
	t52 = (-0.2e1 * t102 / t71 ^ 2 * (qJD(2) * t76 * t88 * t91 + t73 * t78 * t87) * t85 + (t73 * t111 - t72 * t86 + (t80 * t111 + (0.2e1 * t88 * t91 ^ 2 + t86) * t78) * qJD(2)) * t68) * t84;
	t1 = [0, t52, 0, 0, 0, 0; 0, 0.2e1 * (-t54 * t113 - t60 * t82) / t59 ^ 2 * (-t75 * t113 - t77 * t114) + (t54 * t114 * t115 + t74 * t60 + (-t82 * t53 - t54 * t75 - (-(qJD(2) * t110 - t52 * t78 - t56 * t73 + (t56 * t110 - t80) * t55) * t66 - (t55 * t56 * t78 + t72 + (t52 * t92 + (-qJD(2) * t56 - t55) * t91) * t90) * t65) * t101) * t61) / t59, 0, 0, 0, 0; 0, (-t75 / t70 + 0.1e1 / t70 ^ 2 * t82 * t74 * t115 * t118) * t84 / t89, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:17
	% EndTime: 2019-10-09 22:01:17
	% DurationCPUTime: 0.52s
	% Computational Cost: add. (682->54), mult. (2271->134), div. (423->14), fcn. (2956->11), ass. (0->64)
	t135 = sin(pkin(11));
	t137 = cos(pkin(11));
	t140 = sin(qJ(2));
	t138 = cos(pkin(6));
	t142 = cos(qJ(2));
	t155 = t138 * t142;
	t123 = t135 * t140 - t137 * t155;
	t156 = t138 * t140;
	t124 = t135 * t142 + t137 * t156;
	t136 = sin(pkin(6));
	t157 = t136 * t140;
	t114 = atan2(-t124, t157);
	t110 = sin(t114);
	t111 = cos(t114);
	t97 = -t110 * t124 + t111 * t157;
	t94 = 0.1e1 / t97;
	t126 = t135 * t155 + t137 * t140;
	t139 = sin(qJ(4));
	t141 = cos(qJ(4));
	t159 = t135 * t136;
	t109 = t126 * t139 + t141 * t159;
	t105 = 0.1e1 / t109;
	t132 = 0.1e1 / t140;
	t106 = 0.1e1 / t109 ^ 2;
	t133 = 0.1e1 / t140 ^ 2;
	t95 = 0.1e1 / t97 ^ 2;
	t108 = -t126 * t141 + t139 * t159;
	t104 = t108 ^ 2;
	t101 = t104 * t106 + 0.1e1;
	t127 = -t135 * t156 + t137 * t142;
	t120 = t127 * qJD(2);
	t103 = t109 * qJD(4) - t120 * t141;
	t162 = t106 * t108;
	t153 = qJD(4) * t108;
	t102 = t120 * t139 - t153;
	t163 = t102 * t105 * t106;
	t166 = 0.1e1 / t101 ^ 2 * (t103 * t162 - t104 * t163);
	t160 = t133 * t142;
	t152 = t124 * t160;
	t149 = t123 * t132 + t152;
	t121 = t124 ^ 2;
	t131 = 0.1e1 / t136 ^ 2;
	t115 = t121 * t131 * t133 + 0.1e1;
	t112 = 0.1e1 / t115;
	t130 = 0.1e1 / t136;
	t161 = t112 * t130;
	t90 = t149 * t161;
	t165 = t124 * t90;
	t164 = t127 * t95;
	t154 = qJD(2) * t142;
	t150 = t105 * t141 + t139 * t162;
	t134 = t132 * t133;
	t122 = t127 ^ 2;
	t119 = t126 * qJD(2);
	t118 = t124 * qJD(2);
	t117 = t123 * qJD(2);
	t99 = 0.1e1 / t101;
	t96 = t94 * t95;
	t93 = t122 * t95 + 0.1e1;
	t89 = (qJD(2) * t152 + t117 * t132) * t161;
	t87 = (t136 * t142 - t165) * t111 + (-t90 * t157 + t123) * t110;
	t86 = (-t124 * t89 + t136 * t154) * t111 + (-t89 * t157 + t117) * t110;
	t85 = (-0.2e1 * t149 * (-t117 * t124 * t133 - t121 * t134 * t154) * t131 / t115 ^ 2 + (-t117 * t160 + t118 * t132 + (-t123 * t160 + (-0.2e1 * t134 * t142 ^ 2 - t132) * t124) * qJD(2)) * t112) * t130;
	t1 = [0, t85, 0, 0, 0, 0; 0, 0.2e1 * (t126 * t94 + t87 * t164) / t93 ^ 2 * (-t122 * t96 * t86 - t119 * t164) + (t87 * t119 * t95 - t120 * t94 + (0.2e1 * t87 * t127 * t96 + t126 * t95) * t86 + (-(t117 * t90 + t123 * t89 - t124 * t85 + (-t89 * t90 - qJD(2)) * t157) * t111 - (t89 * t165 + t118 + (-t140 * t85 + (-qJD(2) * t90 - t89) * t142) * t136) * t110) * t164) / t93, 0, 0, 0, 0; 0, t150 * t99 * t119 + (0.2e1 * t150 * t166 + ((qJD(4) * t105 + 0.2e1 * t108 * t163) * t139 + (-t103 * t139 + (t102 - t153) * t141) * t106) * t99) * t127, 0, -0.2e1 * t166 + 0.2e1 * (t103 * t106 * t99 + (-t106 * t166 - t99 * t163) * t108) * t108, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:17
	% EndTime: 2019-10-09 22:01:17
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (1083->57), mult. (2577->135), div. (441->14), fcn. (3313->11), ass. (0->71)
	t169 = sin(pkin(11));
	t171 = cos(pkin(11));
	t174 = cos(qJ(2));
	t172 = cos(pkin(6));
	t173 = sin(qJ(2));
	t190 = t172 * t173;
	t154 = t169 * t174 + t171 * t190;
	t170 = sin(pkin(6));
	t191 = t170 * t173;
	t144 = atan2(-t154, t191);
	t140 = sin(t144);
	t141 = cos(t144);
	t127 = -t140 * t154 + t141 * t191;
	t124 = 0.1e1 / t127;
	t189 = t172 * t174;
	t156 = t169 * t189 + t171 * t173;
	t168 = qJ(4) + qJ(5);
	t160 = sin(t168);
	t161 = cos(t168);
	t192 = t169 * t170;
	t139 = t156 * t160 + t161 * t192;
	t135 = 0.1e1 / t139;
	t165 = 0.1e1 / t173;
	t125 = 0.1e1 / t127 ^ 2;
	t136 = 0.1e1 / t139 ^ 2;
	t166 = 0.1e1 / t173 ^ 2;
	t157 = -t169 * t190 + t171 * t174;
	t201 = 0.2e1 * t157;
	t185 = t171 * t189;
	t188 = qJD(2) * t173;
	t147 = -qJD(2) * t185 + t169 * t188;
	t193 = t166 * t174;
	t186 = t154 * t193;
	t151 = t154 ^ 2;
	t163 = 0.1e1 / t170 ^ 2;
	t145 = t151 * t163 * t166 + 0.1e1;
	t142 = 0.1e1 / t145;
	t162 = 0.1e1 / t170;
	t195 = t142 * t162;
	t119 = (qJD(2) * t186 + t147 * t165) * t195;
	t182 = -t140 * t191 - t141 * t154;
	t187 = t141 * t170 * t174;
	t116 = qJD(2) * t187 + t182 * t119 + t140 * t147;
	t200 = t116 * t124 * t125;
	t138 = -t156 * t161 + t160 * t192;
	t134 = t138 ^ 2;
	t131 = t134 * t136 + 0.1e1;
	t150 = t157 * qJD(2);
	t164 = qJD(4) + qJD(5);
	t184 = t164 * t192 - t150;
	t194 = t156 * t164;
	t133 = t160 * t194 + t184 * t161;
	t196 = t136 * t138;
	t132 = -t184 * t160 + t161 * t194;
	t197 = t132 * t135 * t136;
	t199 = (t133 * t196 - t134 * t197) / t131 ^ 2;
	t198 = t125 * t157;
	t183 = t135 * t161 + t160 * t196;
	t153 = t169 * t173 - t185;
	t181 = t153 * t165 + t186;
	t167 = t165 * t166;
	t152 = t157 ^ 2;
	t149 = t156 * qJD(2);
	t148 = t154 * qJD(2);
	t128 = 0.1e1 / t131;
	t123 = t152 * t125 + 0.1e1;
	t120 = t181 * t195;
	t117 = t182 * t120 + t140 * t153 + t187;
	t115 = (-0.2e1 * t181 / t145 ^ 2 * (-qJD(2) * t151 * t167 * t174 - t147 * t154 * t166) * t163 + (-t147 * t193 + t148 * t165 + (-t153 * t193 + (-0.2e1 * t167 * t174 ^ 2 - t165) * t154) * qJD(2)) * t142) * t162;
	t113 = -0.2e1 * t199 + 0.2e1 * (t128 * t133 * t136 + (-t128 * t197 - t136 * t199) * t138) * t138;
	t1 = [0, t115, 0, 0, 0, 0; 0, 0.2e1 * (t117 * t198 + t124 * t156) / t123 ^ 2 * (-t149 * t198 - t152 * t200) + (t117 * t200 * t201 - t150 * t124 + (t156 * t116 + t117 * t149 + (-(-t170 * t188 - t115 * t154 + t120 * t147 + (-t120 * t191 + t153) * t119) * t141 - (t119 * t120 * t154 + t148 + (-t115 * t173 + (-qJD(2) * t120 - t119) * t174) * t170) * t140) * t157) * t125) / t123, 0, 0, 0, 0; 0, t183 * t199 * t201 + (t183 * t149 + ((t135 * t164 + 0.2e1 * t138 * t197) * t160 + (-t133 * t160 + (-t138 * t164 + t132) * t161) * t136) * t157) * t128, 0, t113, t113, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:17
	% EndTime: 2019-10-09 22:01:18
	% DurationCPUTime: 1.66s
	% Computational Cost: add. (9025->112), mult. (13312->233), div. (822->12), fcn. (17117->13), ass. (0->113)
	t252 = qJ(4) + qJ(5);
	t250 = cos(t252);
	t251 = qJD(4) + qJD(5);
	t249 = sin(t252);
	t253 = sin(pkin(11));
	t255 = cos(pkin(11));
	t258 = sin(qJ(2));
	t256 = cos(pkin(6));
	t260 = cos(qJ(2));
	t288 = t256 * t260;
	t272 = -t253 * t258 + t255 * t288;
	t269 = t272 * t249;
	t289 = t256 * t258;
	t245 = t253 * t260 + t255 * t289;
	t254 = sin(pkin(6));
	t292 = t254 * t255;
	t277 = t245 * qJD(2) + t251 * t292;
	t208 = t277 * t250 + t251 * t269;
	t268 = t272 * t250;
	t229 = t249 * t292 - t268;
	t226 = t229 ^ 2;
	t290 = t254 * t260;
	t278 = t250 * t290;
	t238 = t256 * t249 + t278;
	t236 = 0.1e1 / t238 ^ 2;
	t219 = t226 * t236 + 0.1e1;
	t217 = 0.1e1 / t219;
	t291 = t254 * t258;
	t271 = qJD(2) * t291 - t251 * t256;
	t279 = t249 * t290;
	t224 = -t271 * t250 - t251 * t279;
	t235 = 0.1e1 / t238;
	t300 = t229 * t236;
	t189 = (t208 * t235 - t224 * t300) * t217;
	t220 = atan2(t229, t238);
	t215 = sin(t220);
	t216 = cos(t220);
	t275 = -t215 * t238 + t216 * t229;
	t185 = t275 * t189 + t215 * t208 + t216 * t224;
	t199 = t215 * t229 + t216 * t238;
	t196 = 0.1e1 / t199;
	t197 = 0.1e1 / t199 ^ 2;
	t314 = t185 * t196 * t197;
	t246 = t253 * t288 + t255 * t258;
	t293 = t253 * t254;
	t227 = -t246 * t250 + t249 * t293;
	t313 = 0.2e1 * t227 * t314;
	t301 = t224 * t235 * t236;
	t312 = (t208 * t300 - t226 * t301) / t219 ^ 2;
	t281 = t229 * t291;
	t270 = t235 * t245 + t236 * t281;
	t311 = t250 * t270;
	t228 = t246 * t249 + t250 * t293;
	t259 = cos(qJ(6));
	t280 = t253 * t289;
	t247 = t255 * t260 - t280;
	t257 = sin(qJ(6));
	t297 = t247 * t257;
	t214 = t228 * t259 + t297;
	t210 = 0.1e1 / t214;
	t211 = 0.1e1 / t214 ^ 2;
	t287 = qJD(2) * t260;
	t244 = -qJD(2) * t280 + t255 * t287;
	t294 = t250 * t251;
	t205 = t246 * t294 + (-t251 * t293 + t244) * t249;
	t243 = t246 * qJD(2);
	t200 = t214 * qJD(6) + t205 * t257 + t243 * t259;
	t296 = t247 * t259;
	t213 = t228 * t257 - t296;
	t209 = t213 ^ 2;
	t204 = t209 * t211 + 0.1e1;
	t305 = t211 * t213;
	t286 = qJD(6) * t213;
	t201 = t205 * t259 - t243 * t257 - t286;
	t308 = t201 * t210 * t211;
	t310 = (t200 * t305 - t209 * t308) / t204 ^ 2;
	t309 = t197 * t227;
	t206 = t228 * t251 - t244 * t250;
	t307 = t206 * t197;
	t306 = t210 * t257;
	t304 = t213 * t259;
	t303 = t215 * t227;
	t302 = t216 * t227;
	t239 = t256 * t250 - t279;
	t299 = t229 * t239;
	t298 = t247 * t250;
	t295 = t249 * t251;
	t225 = t227 ^ 2;
	t195 = t225 * t197 + 0.1e1;
	t285 = 0.2e1 * (-t225 * t314 + t227 * t307) / t195 ^ 2;
	t284 = -0.2e1 * t310;
	t282 = t213 * t308;
	t276 = qJD(6) * t247 * t249 + t244;
	t274 = t211 * t304 - t306;
	t230 = t250 * t292 + t269;
	t273 = -t230 * t235 + t236 * t299;
	t267 = -qJD(6) * t246 - t243 * t249 + t247 * t294;
	t241 = t272 * qJD(2);
	t223 = t271 * t249 - t251 * t278;
	t222 = -t246 * t257 + t249 * t296;
	t221 = t246 * t259 + t249 * t297;
	t207 = t277 * t249 - t251 * t268;
	t202 = 0.1e1 / t204;
	t193 = 0.1e1 / t195;
	t191 = t217 * t311;
	t190 = t273 * t217;
	t187 = (t215 * t245 - t216 * t291) * t250 + t275 * t191;
	t186 = -t275 * t190 + t215 * t230 + t216 * t239;
	t184 = 0.2e1 * t273 * t312 + (0.2e1 * t299 * t301 - t207 * t235 + (-t208 * t239 - t223 * t229 - t224 * t230) * t236) * t217;
	t182 = -0.2e1 * t311 * t312 + (-t270 * t295 + (-0.2e1 * t281 * t301 + t235 * t241 + (-t224 * t245 + (t208 * t258 + t229 * t287) * t254) * t236) * t250) * t217;
	t181 = t274 * t227 * t284 + (t274 * t206 + ((-qJD(6) * t210 - 0.2e1 * t282) * t259 + (t200 * t259 + (t201 - t286) * t257) * t211) * t227) * t202;
	t180 = (t186 * t309 - t196 * t228) * t285 + (t186 * t313 + t205 * t196 + (-t228 * t185 - t186 * t206 - (t184 * t229 - t190 * t208 + t223 + (t190 * t238 + t230) * t189) * t302 - (-t184 * t238 + t190 * t224 - t207 + (t190 * t229 - t239) * t189) * t303) * t197) * t193;
	t1 = [0, t182, 0, t184, t184, 0; 0, (t187 * t309 + t196 * t298) * t285 + ((t243 * t250 + t247 * t295) * t196 + (-t307 + t313) * t187 + (t298 * t185 - (t182 * t229 + t191 * t208 + (-t250 * t287 + t258 * t295) * t254 + (-t191 * t238 + t245 * t250) * t189) * t302 - (-t245 * t295 - t182 * t238 - t191 * t224 + t241 * t250 + (-t191 * t229 + t250 * t291) * t189) * t303) * t197) * t193, 0, t180, t180, 0; 0, 0.2e1 * (-t210 * t221 + t222 * t305) * t310 + (0.2e1 * t222 * t282 + t276 * t210 * t259 + t267 * t306 + (t276 * t213 * t257 - t222 * t200 - t221 * t201 - t267 * t304) * t211) * t202, 0, t181, t181, t284 + 0.2e1 * (t200 * t211 * t202 + (-t202 * t308 - t211 * t310) * t213) * t213;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end