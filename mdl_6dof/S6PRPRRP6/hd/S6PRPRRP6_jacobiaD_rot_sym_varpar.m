% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRP6
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
%   Wie in S6PRPRRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:52
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRP6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:56
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
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:57
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (536->40), mult. (1574->99), div. (384->14), fcn. (2127->9), ass. (0->50)
	t90 = sin(pkin(6));
	t85 = 0.1e1 / t90 ^ 2;
	t89 = sin(pkin(10));
	t118 = 0.1e1 / t89 ^ 2 * t85;
	t92 = cos(qJ(2));
	t110 = t90 * t92;
	t108 = cos(pkin(10));
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
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:57
	% DurationCPUTime: 0.51s
	% Computational Cost: add. (682->54), mult. (2271->134), div. (423->14), fcn. (2956->11), ass. (0->64)
	t135 = sin(pkin(10));
	t137 = cos(pkin(10));
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
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:58
	% DurationCPUTime: 1.32s
	% Computational Cost: add. (3002->108), mult. (9085->231), div. (559->12), fcn. (11668->13), ass. (0->105)
	t208 = sin(pkin(10));
	t210 = cos(pkin(10));
	t217 = cos(qJ(2));
	t211 = cos(pkin(6));
	t214 = sin(qJ(2));
	t243 = t211 * t214;
	t202 = t208 * t217 + t210 * t243;
	t196 = t202 * qJD(2);
	t209 = sin(pkin(6));
	t216 = cos(qJ(4));
	t213 = sin(qJ(4));
	t242 = t211 * t217;
	t227 = -t208 * t214 + t210 * t242;
	t225 = t227 * t213;
	t167 = qJD(4) * t225 + (t210 * t209 * qJD(4) + t196) * t216;
	t247 = t209 * t213;
	t186 = t210 * t247 - t227 * t216;
	t183 = t186 ^ 2;
	t244 = t209 * t217;
	t205 = t211 * t213 + t216 * t244;
	t200 = 0.1e1 / t205 ^ 2;
	t178 = t183 * t200 + 0.1e1;
	t176 = 0.1e1 / t178;
	t206 = t211 * t216 - t213 * t244;
	t246 = t209 * t214;
	t232 = qJD(2) * t246;
	t189 = t206 * qJD(4) - t216 * t232;
	t199 = 0.1e1 / t205;
	t253 = t186 * t200;
	t148 = (t167 * t199 - t189 * t253) * t176;
	t179 = atan2(t186, t205);
	t174 = sin(t179);
	t175 = cos(t179);
	t230 = -t174 * t205 + t175 * t186;
	t144 = t230 * t148 + t174 * t167 + t175 * t189;
	t158 = t174 * t186 + t175 * t205;
	t155 = 0.1e1 / t158;
	t156 = 0.1e1 / t158 ^ 2;
	t266 = t144 * t155 * t156;
	t203 = t208 * t242 + t210 * t214;
	t184 = -t203 * t216 + t208 * t247;
	t265 = 0.2e1 * t184 * t266;
	t251 = t189 * t199 * t200;
	t264 = (t167 * t253 - t183 * t251) / t178 ^ 2;
	t234 = t186 * t246;
	t226 = t199 * t202 + t200 * t234;
	t263 = t216 * t226;
	t245 = t209 * t216;
	t185 = t203 * t213 + t208 * t245;
	t233 = t208 * t243;
	t204 = t210 * t217 - t233;
	t212 = sin(qJ(5));
	t215 = cos(qJ(5));
	t173 = t185 * t215 + t204 * t212;
	t169 = 0.1e1 / t173;
	t170 = 0.1e1 / t173 ^ 2;
	t241 = qJD(2) * t217;
	t198 = -qJD(2) * t233 + t210 * t241;
	t164 = -t184 * qJD(4) + t198 * t213;
	t197 = t203 * qJD(2);
	t159 = t173 * qJD(5) + t164 * t212 + t197 * t215;
	t249 = t204 * t215;
	t172 = t185 * t212 - t249;
	t168 = t172 ^ 2;
	t163 = t168 * t170 + 0.1e1;
	t257 = t170 * t172;
	t239 = qJD(5) * t172;
	t160 = t164 * t215 - t197 * t212 - t239;
	t260 = t160 * t169 * t170;
	t262 = (t159 * t257 - t168 * t260) / t163 ^ 2;
	t261 = t156 * t184;
	t165 = t185 * qJD(4) - t198 * t216;
	t259 = t165 * t156;
	t258 = t169 * t212;
	t256 = t172 * t215;
	t255 = t174 * t184;
	t254 = t175 * t184;
	t252 = t186 * t206;
	t250 = t204 * t213;
	t248 = t204 * t216;
	t240 = qJD(4) * t213;
	t182 = t184 ^ 2;
	t154 = t182 * t156 + 0.1e1;
	t238 = 0.2e1 * (-t182 * t266 + t184 * t259) / t154 ^ 2;
	t237 = -0.2e1 * t262;
	t235 = t172 * t260;
	t231 = qJD(5) * t250 + t198;
	t229 = t170 * t256 - t258;
	t187 = t210 * t245 + t225;
	t228 = -t187 * t199 + t200 * t252;
	t224 = qJD(4) * t248 - qJD(5) * t203 - t197 * t213;
	t195 = t227 * qJD(2);
	t188 = -t205 * qJD(4) + t213 * t232;
	t181 = -t203 * t212 + t213 * t249;
	t180 = t203 * t215 + t212 * t250;
	t166 = t186 * qJD(4) + t196 * t213;
	t161 = 0.1e1 / t163;
	t151 = 0.1e1 / t154;
	t150 = t176 * t263;
	t149 = t228 * t176;
	t146 = (t174 * t202 - t175 * t246) * t216 + t230 * t150;
	t145 = -t230 * t149 + t174 * t187 + t175 * t206;
	t143 = 0.2e1 * t228 * t264 + (0.2e1 * t251 * t252 - t166 * t199 + (-t167 * t206 - t186 * t188 - t187 * t189) * t200) * t176;
	t141 = -0.2e1 * t263 * t264 + (-t226 * t240 + (-0.2e1 * t234 * t251 + t195 * t199 + (-t189 * t202 + (t167 * t214 + t186 * t241) * t209) * t200) * t216) * t176;
	t1 = [0, t141, 0, t143, 0, 0; 0, (t146 * t261 + t155 * t248) * t238 + ((t197 * t216 + t204 * t240) * t155 + (-t259 + t265) * t146 + (t248 * t144 - (t141 * t186 + t150 * t167 + (t214 * t240 - t216 * t241) * t209 + (-t150 * t205 + t202 * t216) * t148) * t254 - (-t202 * t240 - t141 * t205 - t150 * t189 + t195 * t216 + (-t150 * t186 + t214 * t245) * t148) * t255) * t156) * t151, 0, (t145 * t261 - t155 * t185) * t238 + (t145 * t265 + t164 * t155 + (-t185 * t144 - t145 * t165 - (t143 * t186 - t149 * t167 + t188 + (t149 * t205 + t187) * t148) * t254 - (-t143 * t205 + t149 * t189 - t166 + (t149 * t186 - t206) * t148) * t255) * t156) * t151, 0, 0; 0, 0.2e1 * (-t169 * t180 + t181 * t257) * t262 + (0.2e1 * t181 * t235 + t231 * t169 * t215 + t224 * t258 + (t231 * t172 * t212 - t181 * t159 - t180 * t160 - t224 * t256) * t170) * t161, 0, t229 * t184 * t237 + (t229 * t165 + ((-qJD(5) * t169 - 0.2e1 * t235) * t215 + (t159 * t215 + (t160 - t239) * t212) * t170) * t184) * t161, t237 + 0.2e1 * (t159 * t170 * t161 + (-t161 * t260 - t170 * t262) * t172) * t172, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:57
	% EndTime: 2019-10-09 21:51:59
	% DurationCPUTime: 2.45s
	% Computational Cost: add. (7242->157), mult. (21125->305), div. (784->12), fcn. (27225->13), ass. (0->124)
	t256 = sin(qJ(4));
	t259 = cos(qJ(4));
	t257 = sin(qJ(2));
	t260 = cos(qJ(2));
	t320 = cos(pkin(10));
	t321 = cos(pkin(6));
	t284 = t321 * t320;
	t319 = sin(pkin(10));
	t273 = -t319 * t257 + t260 * t284;
	t254 = sin(pkin(6));
	t290 = t254 * t320;
	t229 = t256 * t290 - t259 * t273;
	t274 = t257 * t284 + t319 * t260;
	t245 = t274 * qJD(2);
	t211 = t229 * qJD(4) + t245 * t256;
	t258 = cos(qJ(5));
	t255 = sin(qJ(5));
	t271 = t274 * t255;
	t275 = t256 * t273 + t259 * t290;
	t222 = -t258 * t275 + t271;
	t269 = t273 * qJD(2);
	t194 = t222 * qJD(5) + t211 * t255 - t258 * t269;
	t270 = t274 * t258;
	t220 = -t255 * t275 - t270;
	t217 = t220 ^ 2;
	t306 = t254 * t260;
	t276 = t256 * t306 - t321 * t259;
	t303 = t257 * t258;
	t277 = t254 * t303 + t255 * t276;
	t232 = 0.1e1 / t277 ^ 2;
	t204 = t217 * t232 + 0.1e1;
	t202 = 0.1e1 / t204;
	t251 = -t321 * t256 - t259 * t306;
	t302 = qJD(2) * t254;
	t292 = t257 * t302;
	t234 = t251 * qJD(4) + t256 * t292;
	t305 = t255 * t257;
	t237 = t254 * t305 - t258 * t276;
	t291 = t260 * t302;
	t206 = t237 * qJD(5) + t234 * t255 - t258 * t291;
	t231 = 0.1e1 / t277;
	t310 = t220 * t232;
	t179 = (t194 * t231 + t206 * t310) * t202;
	t205 = atan2(-t220, -t277);
	t199 = sin(t205);
	t200 = cos(t205);
	t282 = t199 * t277 - t200 * t220;
	t175 = t282 * t179 - t199 * t194 + t200 * t206;
	t191 = -t199 * t220 - t200 * t277;
	t188 = 0.1e1 / t191;
	t189 = 0.1e1 / t191 ^ 2;
	t326 = t175 * t188 * t189;
	t283 = t321 * t319;
	t272 = -t320 * t257 - t260 * t283;
	t289 = t254 * t319;
	t228 = -t256 * t272 + t259 * t289;
	t250 = -t257 * t283 + t320 * t260;
	t218 = t228 * t255 - t250 * t258;
	t325 = -0.2e1 * t218;
	t287 = 0.2e1 * t218 * t326;
	t278 = t229 * t231 + t251 * t310;
	t324 = t255 * t278;
	t323 = qJD(5) * t270 + t255 * t269;
	t312 = t206 * t231 * t232;
	t322 = -0.2e1 * (t194 * t310 + t217 * t312) / t204 ^ 2;
	t219 = t228 * t258 + t250 * t255;
	t214 = 0.1e1 / t219;
	t215 = 0.1e1 / t219 ^ 2;
	t285 = t256 * t289;
	t227 = -t259 * t272 - t285;
	t226 = t227 ^ 2;
	t309 = t226 * t215;
	t201 = 0.1e1 + t309;
	t247 = t250 * qJD(2);
	t210 = -t228 * qJD(4) + t247 * t259;
	t301 = qJD(4) * t259;
	t209 = qJD(4) * t285 - t247 * t256 + t272 * t301;
	t246 = t272 * qJD(2);
	t193 = -t218 * qJD(5) - t209 * t258 + t246 * t255;
	t315 = t193 * t214 * t215;
	t293 = t226 * t315;
	t311 = t215 * t227;
	t318 = (t210 * t311 - t293) / t201 ^ 2;
	t317 = t189 * t218;
	t192 = t219 * qJD(5) - t209 * t255 - t246 * t258;
	t316 = t192 * t189;
	t314 = t199 * t218;
	t313 = t200 * t218;
	t308 = t227 * t255;
	t307 = t250 * t256;
	t304 = t256 * t258;
	t300 = qJD(5) * t255;
	t299 = qJD(5) * t256;
	t298 = qJD(5) * t258;
	t213 = t218 ^ 2;
	t187 = t213 * t189 + 0.1e1;
	t297 = 0.2e1 * (-t213 * t326 + t218 * t316) / t187 ^ 2;
	t296 = 0.2e1 * t318;
	t294 = t227 * t315;
	t286 = 0.2e1 * t220 * t312;
	t280 = t222 * t231 + t237 * t310;
	t223 = t256 * t271 - t258 * t273;
	t240 = (t256 * t305 - t258 * t260) * t254;
	t279 = t223 * t231 + t240 * t310;
	t235 = t276 * qJD(4) + t259 * t292;
	t225 = t250 * t304 + t255 * t272;
	t224 = t255 * t307 - t258 * t272;
	t212 = t275 * qJD(4) + t245 * t259;
	t208 = ((qJD(2) + t299) * t303 + (t257 * t301 + (qJD(2) * t256 + qJD(5)) * t260) * t255) * t254;
	t207 = t277 * qJD(5) + t234 * t258 + t255 * t291;
	t197 = 0.1e1 / t201;
	t196 = t245 * t258 + t323 * t256 + t271 * t301 + t273 * t300;
	t195 = t211 * t258 + t275 * t300 + t323;
	t185 = 0.1e1 / t187;
	t184 = t202 * t324;
	t183 = t279 * t202;
	t181 = t280 * t202;
	t178 = (-t199 * t229 + t200 * t251) * t255 + t282 * t184;
	t177 = t282 * t183 - t199 * t223 + t200 * t240;
	t176 = t282 * t181 - t199 * t222 + t200 * t237;
	t174 = t279 * t322 + (t240 * t286 + t196 * t231 + (t194 * t240 + t206 * t223 + t208 * t220) * t232) * t202;
	t172 = t280 * t322 + (t237 * t286 + t195 * t231 + (t194 * t237 + t206 * t222 + t207 * t220) * t232) * t202;
	t171 = t322 * t324 + (t278 * t298 + (t251 * t286 + t212 * t231 + (t194 * t251 + t206 * t229 + t220 * t235) * t232) * t255) * t202;
	t1 = [0, t174, 0, t171, t172, 0; 0, (t177 * t317 - t188 * t224) * t297 + (t177 * t287 + (-t224 * t175 - t177 * t192 - (-t174 * t220 - t183 * t194 + t208 + (t183 * t277 - t223) * t179) * t313 - (t174 * t277 - t183 * t206 - t196 + (t183 * t220 - t240) * t179) * t314) * t189 + ((t250 * t299 + t247) * t258 + (qJD(5) * t272 + t246 * t256 + t250 * t301) * t255) * t188) * t185, 0, (t178 * t317 - t188 * t308) * t297 + ((t210 * t255 + t227 * t298) * t188 + (-t316 + t287) * t178 + (-t308 * t175 - (t251 * t298 - t171 * t220 - t184 * t194 + t235 * t255 + (t184 * t277 - t229 * t255) * t179) * t313 - (-t229 * t298 + t171 * t277 - t184 * t206 - t212 * t255 + (t184 * t220 - t251 * t255) * t179) * t314) * t189) * t185, (t176 * t317 - t188 * t219) * t297 + (t176 * t287 + t193 * t188 + (-t219 * t175 - t176 * t192 - (-t172 * t220 - t181 * t194 + t207 + (t181 * t277 - t222) * t179) * t313 - (t172 * t277 - t181 * t206 - t195 + (t181 * t220 - t237) * t179) * t314) * t189) * t185, 0; 0, (-t214 * t250 * t259 + t225 * t311) * t296 + (0.2e1 * t225 * t294 + (-qJD(4) * t307 + t246 * t259) * t214 + (-(t246 * t304 - t247 * t255 + t272 * t298) * t227 - t225 * t210 + (-t259 * t193 - (-t255 * t299 + t258 * t301) * t227) * t250) * t215) * t197, 0, (t214 * t228 + t258 * t309) * t296 + (0.2e1 * t258 * t293 + t209 * t214 + (-0.2e1 * t210 * t227 * t258 + t193 * t228 + t226 * t300) * t215) * t197, t311 * t318 * t325 + (t294 * t325 + (t192 * t227 + t210 * t218) * t215) * t197, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end