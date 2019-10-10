% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR7
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
%   Wie in S6PRPRPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:40
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
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
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
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
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:57
	% DurationCPUTime: 0.52s
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
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:57
	% DurationCPUTime: 1.01s
	% Computational Cost: add. (2385->85), mult. (7313->189), div. (516->12), fcn. (9465->11), ass. (0->87)
	t158 = sin(pkin(10));
	t160 = cos(pkin(10));
	t165 = cos(qJ(2));
	t161 = cos(pkin(6));
	t163 = sin(qJ(2));
	t185 = t161 * t163;
	t153 = t158 * t165 + t160 * t185;
	t144 = t153 * qJD(2);
	t159 = sin(pkin(6));
	t164 = cos(qJ(4));
	t162 = sin(qJ(4));
	t184 = t161 * t165;
	t174 = -t158 * t163 + t160 * t184;
	t172 = t174 * t162;
	t121 = qJD(4) * t172 + (qJD(4) * t159 * t160 + t144) * t164;
	t189 = t159 * t162;
	t136 = t160 * t189 - t164 * t174;
	t133 = t136 ^ 2;
	t186 = t159 * t165;
	t156 = t161 * t162 + t164 * t186;
	t151 = 0.1e1 / t156 ^ 2;
	t129 = t133 * t151 + 0.1e1;
	t126 = 0.1e1 / t129;
	t157 = t161 * t164 - t162 * t186;
	t188 = t159 * t163;
	t177 = qJD(2) * t188;
	t139 = qJD(4) * t157 - t164 * t177;
	t150 = 0.1e1 / t156;
	t195 = t136 * t151;
	t106 = (t121 * t150 - t139 * t195) * t126;
	t130 = atan2(t136, t156);
	t122 = sin(t130);
	t123 = cos(t130);
	t176 = -t122 * t156 + t123 * t136;
	t103 = t106 * t176 + t122 * t121 + t123 * t139;
	t117 = t122 * t136 + t123 * t156;
	t114 = 0.1e1 / t117;
	t115 = 0.1e1 / t117 ^ 2;
	t204 = t103 * t114 * t115;
	t154 = t158 * t184 + t160 * t163;
	t134 = -t154 * t164 + t158 * t189;
	t203 = 0.2e1 * t134 * t204;
	t193 = t139 * t150 * t151;
	t202 = (t121 * t195 - t133 * t193) / t129 ^ 2;
	t145 = t154 * qJD(2);
	t155 = -t158 * t185 + t160 * t165;
	t148 = 0.1e1 / t155 ^ 2;
	t201 = t145 * t148;
	t178 = t136 * t188;
	t173 = t150 * t153 + t151 * t178;
	t200 = t164 * t173;
	t147 = 0.1e1 / t155;
	t199 = t115 * t134;
	t187 = t159 * t164;
	t135 = t154 * t162 + t158 * t187;
	t146 = t155 * qJD(2);
	t119 = qJD(4) * t135 - t146 * t164;
	t198 = t119 * t115;
	t197 = t122 * t134;
	t196 = t123 * t134;
	t194 = t136 * t157;
	t192 = t147 * t201;
	t191 = t148 * t154;
	t190 = t155 * t164;
	t183 = qJD(2) * t165;
	t182 = qJD(4) * t162;
	t131 = t134 ^ 2;
	t112 = t115 * t131 + 0.1e1;
	t181 = 0.2e1 * (-t131 * t204 + t134 * t198) / t112 ^ 2;
	t118 = -qJD(4) * t134 + t146 * t162;
	t132 = t135 ^ 2;
	t128 = t132 * t148 + 0.1e1;
	t180 = 0.2e1 * (t118 * t135 * t148 + t132 * t192) / t128 ^ 2;
	t137 = t160 * t187 + t172;
	t175 = -t137 * t150 + t151 * t194;
	t143 = t174 * qJD(2);
	t138 = -qJD(4) * t156 + t162 * t177;
	t124 = 0.1e1 / t128;
	t120 = qJD(4) * t136 + t144 * t162;
	t109 = 0.1e1 / t112;
	t108 = t126 * t200;
	t107 = t175 * t126;
	t105 = (t122 * t153 - t123 * t188) * t164 + t176 * t108;
	t104 = -t107 * t176 + t122 * t137 + t123 * t157;
	t102 = 0.2e1 * t175 * t202 + (0.2e1 * t193 * t194 - t120 * t150 + (-t121 * t157 - t136 * t138 - t137 * t139) * t151) * t126;
	t100 = -0.2e1 * t200 * t202 + (-t173 * t182 + (-0.2e1 * t178 * t193 + t143 * t150 + (-t139 * t153 + (t121 * t163 + t136 * t183) * t159) * t151) * t164) * t126;
	t1 = [0, t100, 0, t102, 0, 0; 0, (t105 * t199 + t114 * t190) * t181 + ((t145 * t164 + t155 * t182) * t114 + (-t198 + t203) * t105 + (t190 * t103 - (t100 * t136 + t108 * t121 + (t163 * t182 - t164 * t183) * t159 + (-t108 * t156 + t153 * t164) * t106) * t196 - (-t153 * t182 - t100 * t156 - t108 * t139 + t143 * t164 + (-t108 * t136 + t163 * t187) * t106) * t197) * t115) * t109, 0, (t104 * t199 - t114 * t135) * t181 + (t104 * t203 + t118 * t114 + (-t135 * t103 - t104 * t119 - (t102 * t136 - t107 * t121 + t138 + (t107 * t156 + t137) * t106) * t196 - (-t102 * t156 + t107 * t139 - t120 + (t107 * t136 - t157) * t106) * t197) * t115) * t109, 0, 0; 0, (-t135 * t191 - t162) * t180 + (t118 * t191 + qJD(4) * t164 + (t146 * t148 + 0.2e1 * t154 * t192) * t135) * t124, 0, t134 * t147 * t180 + (-t119 * t147 - t134 * t201) * t124, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:58
	% DurationCPUTime: 1.38s
	% Computational Cost: add. (3002->108), mult. (9085->223), div. (559->12), fcn. (11668->13), ass. (0->102)
	t207 = cos(pkin(10));
	t213 = cos(qJ(4));
	t205 = sin(pkin(10));
	t211 = sin(qJ(2));
	t208 = cos(pkin(6));
	t214 = cos(qJ(2));
	t240 = t208 * t214;
	t224 = -t205 * t211 + t207 * t240;
	t206 = sin(pkin(6));
	t210 = sin(qJ(4));
	t245 = t206 * t210;
	t186 = t207 * t245 - t224 * t213;
	t241 = t208 * t211;
	t199 = t205 * t214 + t207 * t241;
	t193 = t199 * qJD(2);
	t165 = qJD(4) * t186 + t193 * t210;
	t243 = t206 * t213;
	t187 = t207 * t243 + t224 * t210;
	t183 = t187 ^ 2;
	t242 = t206 * t214;
	t223 = -t208 * t213 + t210 * t242;
	t197 = 0.1e1 / t223 ^ 2;
	t177 = t183 * t197 + 0.1e1;
	t175 = 0.1e1 / t177;
	t202 = -t208 * t210 - t213 * t242;
	t244 = t206 * t211;
	t232 = qJD(2) * t244;
	t188 = qJD(4) * t202 + t210 * t232;
	t196 = 0.1e1 / t223;
	t249 = t187 * t197;
	t147 = (t165 * t196 - t188 * t249) * t175;
	t178 = atan2(t187, -t223);
	t173 = sin(t178);
	t174 = cos(t178);
	t228 = t173 * t223 + t174 * t187;
	t143 = t228 * t147 - t173 * t165 + t174 * t188;
	t157 = t173 * t187 - t174 * t223;
	t154 = 0.1e1 / t157;
	t155 = 0.1e1 / t157 ^ 2;
	t262 = t143 * t154 * t155;
	t200 = t205 * t240 + t207 * t211;
	t184 = -t200 * t213 + t205 * t245;
	t233 = t205 * t241;
	t201 = t207 * t214 - t233;
	t209 = sin(qJ(6));
	t212 = cos(qJ(6));
	t227 = t184 * t212 - t201 * t209;
	t261 = t227 * qJD(6);
	t185 = t200 * t210 + t205 * t243;
	t260 = 0.2e1 * t185 * t262;
	t222 = -t196 * t199 + t244 * t249;
	t259 = t210 * t222;
	t172 = t184 * t209 + t201 * t212;
	t168 = 0.1e1 / t172;
	t169 = 0.1e1 / t172 ^ 2;
	t239 = qJD(2) * t214;
	t195 = -qJD(2) * t233 + t207 * t239;
	t164 = t185 * qJD(4) - t195 * t213;
	t194 = t200 * qJD(2);
	t158 = t172 * qJD(6) - t164 * t212 - t194 * t209;
	t167 = t227 ^ 2;
	t162 = t167 * t169 + 0.1e1;
	t253 = t169 * t227;
	t159 = t164 * t209 - t194 * t212 + t261;
	t256 = t159 * t168 * t169;
	t258 = (-t158 * t253 - t167 * t256) / t162 ^ 2;
	t257 = t155 * t185;
	t163 = -t184 * qJD(4) + t195 * t210;
	t255 = t163 * t155;
	t254 = t168 * t212;
	t252 = t227 * t209;
	t251 = t173 * t185;
	t250 = t174 * t185;
	t248 = t188 * t196 * t197;
	t247 = t201 * t210;
	t246 = t201 * t213;
	t238 = qJD(4) * t213;
	t182 = t185 ^ 2;
	t153 = t182 * t155 + 0.1e1;
	t237 = 0.2e1 * (-t182 * t262 + t185 * t255) / t153 ^ 2;
	t236 = 0.2e1 * t258;
	t235 = 0.2e1 * (-t165 * t249 + t183 * t248) / t177 ^ 2;
	t231 = -0.2e1 * t227 * t256;
	t230 = -0.2e1 * t187 * t248;
	t229 = -qJD(6) * t246 - t195;
	t226 = -t169 * t252 + t254;
	t225 = -t186 * t196 + t202 * t249;
	t221 = qJD(4) * t247 + qJD(6) * t200 + t194 * t213;
	t192 = t224 * qJD(2);
	t189 = qJD(4) * t223 + t213 * t232;
	t180 = -t200 * t212 - t209 * t246;
	t179 = -t200 * t209 + t212 * t246;
	t166 = qJD(4) * t187 + t193 * t213;
	t160 = 0.1e1 / t162;
	t150 = 0.1e1 / t153;
	t149 = t175 * t259;
	t148 = t225 * t175;
	t145 = (-t173 * t199 + t174 * t244) * t210 - t228 * t149;
	t144 = -t228 * t148 - t173 * t186 + t174 * t202;
	t142 = t225 * t235 + (t202 * t230 + t166 * t196 + (t165 * t202 + t186 * t188 - t187 * t189) * t197) * t175;
	t140 = t235 * t259 + (-t222 * t238 + (t230 * t244 + t192 * t196 + (t188 * t199 + (t165 * t211 - t187 * t239) * t206) * t197) * t210) * t175;
	t1 = [0, t140, 0, t142, 0, 0; 0, (t145 * t257 - t154 * t247) * t237 + ((-t194 * t210 + t201 * t238) * t154 + (-t255 + t260) * t145 + (-t247 * t143 - (t140 * t187 + t149 * t165 + (t210 * t239 + t211 * t238) * t206 + (-t149 * t223 - t199 * t210) * t147) * t250 - (-t199 * t238 + t140 * t223 + t149 * t188 - t192 * t210 + (t149 * t187 - t210 * t244) * t147) * t251) * t155) * t150, 0, (t144 * t257 + t154 * t184) * t237 + (t144 * t260 - t164 * t154 + (t184 * t143 - t144 * t163 - (t142 * t187 + t148 * t165 + t189 + (-t148 * t223 - t186) * t147) * t250 - (t142 * t223 + t148 * t188 - t166 + (t148 * t187 - t202) * t147) * t251) * t155) * t150, 0, 0; 0, (-t168 * t179 - t180 * t253) * t236 + (t180 * t231 + t229 * t168 * t209 - t221 * t254 + (t212 * t227 * t229 - t180 * t158 - t179 * t159 + t221 * t252) * t169) * t160, 0, t226 * t185 * t236 + (-t226 * t163 + ((qJD(6) * t168 + t231) * t209 + (-t158 * t209 + (t159 + t261) * t212) * t169) * t185) * t160, 0, -0.2e1 * t258 - 0.2e1 * (t158 * t169 * t160 - (-t160 * t256 - t169 * t258) * t227) * t227;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end