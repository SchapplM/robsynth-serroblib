% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR12
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
%   Wie in S6RPRPRR12_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:06
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR12_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR12_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:05:58
	% EndTime: 2019-10-10 01:05:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:05:58
	% EndTime: 2019-10-10 01:05:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:05:58
	% EndTime: 2019-10-10 01:05:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:05:58
	% EndTime: 2019-10-10 01:05:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:05:58
	% EndTime: 2019-10-10 01:05:59
	% DurationCPUTime: 0.70s
	% Computational Cost: add. (585->68), mult. (1839->158), div. (436->14), fcn. (2165->7), ass. (0->74)
	t85 = sin(qJ(1));
	t116 = qJD(1) * t85;
	t87 = cos(qJ(1));
	t137 = 0.2e1 * t87;
	t77 = t85 ^ 2;
	t82 = 0.1e1 / t87 ^ 2;
	t119 = t77 * t82;
	t84 = sin(qJ(3));
	t72 = t84 ^ 2;
	t70 = t72 * t119 + 0.1e1;
	t80 = t87 ^ 2;
	t81 = 0.1e1 / t87;
	t95 = (t77 / t80 + 0.1e1) * t81 * t116;
	t113 = qJD(3) * t84;
	t86 = cos(qJ(3));
	t98 = t77 * t86 * t113;
	t136 = -0.2e1 * (t72 * t95 + t82 * t98) / t70 ^ 2;
	t111 = qJD(3) * t87;
	t115 = qJD(1) * t86;
	t105 = t85 * t115;
	t74 = 0.1e1 / t84 ^ 2;
	t79 = t86 ^ 2;
	t121 = t74 * t79;
	t71 = t80 * t121 + 0.1e1;
	t68 = 0.1e1 / t71;
	t73 = 0.1e1 / t84;
	t52 = ((t111 * t84 + t105) * t73 + t111 * t121) * t68;
	t134 = -t52 + t111;
	t117 = t87 * t86;
	t65 = atan2(-t117, t84);
	t63 = sin(t65);
	t64 = cos(t65);
	t59 = -t63 * t117 + t64 * t84;
	t56 = 0.1e1 / t59;
	t57 = 0.1e1 / t59 ^ 2;
	t133 = t68 - 0.1e1;
	t120 = t77 * t79;
	t124 = t64 * t86;
	t48 = (-t52 * t87 + qJD(3)) * t124 + (t134 * t84 + t105) * t63;
	t131 = t48 * t56 * t57;
	t55 = t57 * t120 + 0.1e1;
	t114 = qJD(1) * t87;
	t99 = t79 * t85 * t114;
	t132 = (-t120 * t131 + (-t98 + t99) * t57) / t55 ^ 2;
	t130 = t52 * t86;
	t129 = t57 * t85;
	t128 = t57 * t86;
	t122 = t73 * t86;
	t78 = t86 * t79;
	t94 = qJD(3) * (-t73 / t72 * t78 - t122);
	t127 = (-t74 * t99 + t80 * t94) / t71 ^ 2;
	t125 = t63 * t87;
	t123 = t73 * t79;
	t118 = t85 * t86;
	t112 = qJD(3) * t85;
	t110 = -0.2e1 * t131;
	t109 = t86 * t132;
	t108 = t57 * t118;
	t107 = t86 * t127;
	t106 = t68 * t123;
	t104 = 0.1e1 + t121;
	t103 = 0.1e1 + t119;
	t102 = t127 * t137;
	t101 = t87 * t106;
	t100 = t133 * t86 * t63;
	t97 = t104 * t85;
	t96 = t103 * t86;
	t66 = 0.1e1 / t70;
	t62 = t104 * t87 * t68;
	t53 = 0.1e1 / t55;
	t51 = (-t101 * t64 - t100) * t85;
	t50 = t84 * t125 + t124 + (-t64 * t117 - t63 * t84) * t62;
	t49 = -t104 * t102 + (-qJD(1) * t97 + t94 * t137) * t68;
	t1 = [-0.2e1 * t85 * t73 * t107 + (-qJD(3) * t97 + t114 * t122) * t68, 0, t49, 0, 0, 0; (0.2e1 * t56 * t109 + (t56 * t113 + (qJD(1) * t51 + t48) * t128) * t53) * t87 + (-0.2e1 * t57 * t109 * t51 + (((t52 * t101 + t133 * t113 + 0.2e1 * t107) * t63 + (t102 * t123 + t130 + (-t130 + (t74 * t78 + 0.2e1 * t86) * t111) * t68) * t64) * t108 + (t110 * t86 - t113 * t57) * t51 + (t56 + ((t77 - t80) * t64 * t106 - t87 * t100) * t57) * t115) * t53) * t85, 0, 0.2e1 * (-t50 * t128 - t56 * t84) * t85 * t132 + ((t56 * t114 + (-qJD(3) * t50 - t48) * t129) * t84 + (t56 * t112 + (-t49 * t64 * t87 + t134 * t63 + (-qJD(3) * t63 + t116 * t64 + t125 * t52) * t62) * t108 + (t110 * t85 + t57 * t114) * t50 + ((-t49 - t116) * t63 + ((t62 * t87 - 0.1e1) * qJD(3) + (-t62 + t87) * t52) * t64) * t84 * t129) * t86) * t53, 0, 0, 0; t103 * t84 * t136 + (qJD(3) * t96 + 0.2e1 * t84 * t95) * t66, 0, t81 * t118 * t136 + (-t112 * t81 * t84 + qJD(1) * t96) * t66, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:05:58
	% EndTime: 2019-10-10 01:05:59
	% DurationCPUTime: 0.97s
	% Computational Cost: add. (624->92), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->93)
	t129 = cos(qJ(1));
	t123 = t129 ^ 2;
	t125 = sin(qJ(3));
	t118 = t125 ^ 2;
	t128 = cos(qJ(3));
	t121 = 0.1e1 / t128 ^ 2;
	t171 = t118 * t121;
	t114 = t123 * t171 + 0.1e1;
	t117 = t125 * t118;
	t120 = 0.1e1 / t128;
	t170 = t120 * t125;
	t138 = qJD(3) * (t117 * t120 * t121 + t170);
	t126 = sin(qJ(1));
	t163 = qJD(1) * t129;
	t153 = t126 * t163;
	t180 = 0.1e1 / t114 ^ 2 * (t123 * t138 - t153 * t171);
	t192 = -0.2e1 * t180;
	t111 = 0.1e1 / t114;
	t148 = 0.1e1 + t171;
	t190 = t129 * t148;
	t99 = t111 * t190;
	t191 = -t129 * t99 + 0.1e1;
	t144 = qJD(1) * t128 + qJD(5);
	t160 = qJD(3) * t129;
	t189 = t125 * t160 + t144 * t126;
	t166 = t129 * t125;
	t113 = atan2(t166, t128);
	t109 = sin(t113);
	t110 = cos(t113);
	t95 = t109 * t166 + t110 * t128;
	t92 = 0.1e1 / t95;
	t124 = sin(qJ(5));
	t127 = cos(qJ(5));
	t165 = t129 * t127;
	t168 = t126 * t128;
	t108 = -t124 * t168 + t165;
	t102 = 0.1e1 / t108;
	t103 = 0.1e1 / t108 ^ 2;
	t93 = 0.1e1 / t95 ^ 2;
	t188 = t111 - 0.1e1;
	t119 = t126 ^ 2;
	t161 = qJD(3) * t128;
	t157 = t93 * t161;
	t164 = qJD(1) * t126;
	t154 = t125 * t164;
	t173 = t110 * t125;
	t152 = t121 * t160;
	t86 = ((t128 * t160 - t154) * t120 + t118 * t152) * t111;
	t81 = (t129 * t86 - qJD(3)) * t173 + (-t154 + (-t86 + t160) * t128) * t109;
	t186 = t81 * t92 * t93;
	t91 = t119 * t118 * t93 + 0.1e1;
	t187 = (t119 * t125 * t157 + (-t119 * t186 + t93 * t153) * t118) / t91 ^ 2;
	t167 = t129 * t124;
	t107 = t127 * t168 + t167;
	t101 = t107 ^ 2;
	t100 = t101 * t103 + 0.1e1;
	t175 = t103 * t107;
	t145 = qJD(5) * t128 + qJD(1);
	t162 = qJD(3) * t126;
	t169 = t126 * t127;
	t88 = -t145 * t169 + (t125 * t162 - t144 * t129) * t124;
	t182 = t102 * t103 * t88;
	t155 = t128 * t165;
	t87 = -qJD(1) * t155 - qJD(5) * t165 + (qJD(3) * t125 * t127 + t145 * t124) * t126;
	t185 = (-t101 * t182 - t87 * t175) / t100 ^ 2;
	t89 = 0.1e1 / t91;
	t184 = t89 * t93;
	t183 = t92 * t89;
	t179 = t126 * t93;
	t177 = qJD(3) * t99;
	t176 = t102 * t127;
	t174 = t107 * t124;
	t172 = t118 * t120;
	t159 = 0.2e1 * t186;
	t158 = 0.2e1 * t185;
	t156 = t129 * t172;
	t150 = -0.2e1 * t92 * t187;
	t149 = 0.2e1 * t93 * t187;
	t147 = 0.2e1 * t107 * t182;
	t146 = 0.2e1 * t125 * t180;
	t143 = t111 * t156;
	t142 = t188 * t125 * t109;
	t141 = t148 * t126;
	t140 = t145 * t129;
	t139 = t103 * t174 + t176;
	t97 = 0.1e1 / t100;
	t137 = t139 * t97;
	t106 = -t128 * t167 - t169;
	t105 = -t126 * t124 + t155;
	t85 = (-t110 * t143 + t142) * t126;
	t84 = -t191 * t173 + (t129 - t99) * t128 * t109;
	t82 = t190 * t192 + (-qJD(1) * t141 + 0.2e1 * t129 * t138) * t111;
	t1 = [t126 * t120 * t146 + (-qJD(3) * t141 - t163 * t170) * t111, 0, t82, 0, 0, 0; (t161 * t183 + (t150 + (-qJD(1) * t85 - t81) * t184) * t125) * t129 + (t85 * t149 * t125 + (-t85 * t157 + (t85 * t159 + ((-t86 * t143 - t188 * t161 + t146) * t109 + (t156 * t192 + t125 * t86 + (t117 * t152 - (t86 - 0.2e1 * t160) * t125) * t111) * t110) * t179) * t125 + (-t92 + (-(t119 - t123) * t111 * t110 * t172 - t129 * t142) * t93) * t125 * qJD(1)) * t89) * t126, 0, (t163 * t183 + (t150 + (-qJD(3) * t84 - t81) * t184) * t126) * t128 + (t84 * t126 * t149 + (-t92 * t162 - ((t129 * t82 - t164 * t99) * t110 + (t191 * t86 - t160 + t177) * t109) * t125 * t179 + (t126 * t159 - t93 * t163) * t84) * t89 - ((-t82 - t164) * t109 + (-t86 * t99 - qJD(3) + (t86 + t177) * t129) * t110) * t168 * t184) * t125, 0, 0, 0; (-t102 * t105 + t106 * t175) * t158 + (t106 * t147 - t102 * t124 * t140 - t189 * t176 + (t107 * t127 * t140 - t105 * t88 + t106 * t87 - t189 * t174) * t103) * t97, 0, -t126 * t137 * t161 + (-t137 * t163 + (t139 * t158 + ((qJD(5) * t102 + t147) * t124 + (t124 * t87 + (-qJD(5) * t107 + t88) * t127) * t103) * t97) * t126) * t125, 0, -0.2e1 * t185 + 0.2e1 * (-t87 * t103 * t97 + (-t103 * t185 - t97 * t182) * t107) * t107, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:05:58
	% EndTime: 2019-10-10 01:05:59
	% DurationCPUTime: 1.02s
	% Computational Cost: add. (1192->94), mult. (2734->207), div. (498->12), fcn. (3199->9), ass. (0->95)
	t155 = sin(qJ(3));
	t148 = t155 ^ 2;
	t157 = cos(qJ(3));
	t151 = 0.1e1 / t157 ^ 2;
	t199 = t148 * t151;
	t158 = cos(qJ(1));
	t219 = 0.2e1 * t158;
	t218 = t155 * t199;
	t153 = t158 ^ 2;
	t143 = t153 * t199 + 0.1e1;
	t140 = 0.1e1 / t143;
	t150 = 0.1e1 / t157;
	t156 = sin(qJ(1));
	t191 = qJD(1) * t156;
	t179 = t155 * t191;
	t187 = qJD(3) * t158;
	t115 = ((t157 * t187 - t179) * t150 + t187 * t199) * t140;
	t217 = t115 - t187;
	t146 = qJD(5) + qJD(6);
	t173 = qJD(1) * t157 + t146;
	t216 = t155 * t187 + t173 * t156;
	t192 = t158 * t155;
	t142 = atan2(t192, t157);
	t137 = sin(t142);
	t138 = cos(t142);
	t125 = t137 * t192 + t138 * t157;
	t122 = 0.1e1 / t125;
	t154 = qJ(5) + qJ(6);
	t144 = sin(t154);
	t145 = cos(t154);
	t193 = t158 * t145;
	t195 = t156 * t157;
	t135 = -t144 * t195 + t193;
	t129 = 0.1e1 / t135;
	t123 = 0.1e1 / t125 ^ 2;
	t130 = 0.1e1 / t135 ^ 2;
	t215 = 0.2e1 * t156;
	t214 = t140 - 0.1e1;
	t149 = t156 ^ 2;
	t201 = t148 * t149;
	t120 = t123 * t201 + 0.1e1;
	t190 = qJD(1) * t158;
	t170 = t148 * t156 * t190;
	t188 = qJD(3) * t157;
	t202 = t138 * t155;
	t109 = (t115 * t158 - qJD(3)) * t202 + (-t217 * t157 - t179) * t137;
	t211 = t109 * t122 * t123;
	t213 = (-t201 * t211 + (t149 * t155 * t188 + t170) * t123) / t120 ^ 2;
	t174 = t146 * t157 + qJD(1);
	t180 = t157 * t193;
	t113 = -qJD(1) * t180 - t146 * t193 + (qJD(3) * t145 * t155 + t174 * t144) * t156;
	t194 = t158 * t144;
	t134 = t145 * t195 + t194;
	t128 = t134 ^ 2;
	t121 = t128 * t130 + 0.1e1;
	t205 = t130 * t134;
	t189 = qJD(3) * t156;
	t196 = t156 * t145;
	t114 = -t174 * t196 + (t155 * t189 - t173 * t158) * t144;
	t210 = t114 * t129 * t130;
	t212 = (-t113 * t205 - t128 * t210) / t121 ^ 2;
	t209 = t115 * t155;
	t208 = t123 * t155;
	t198 = t150 * t155;
	t166 = qJD(3) * (t150 * t218 + t198);
	t207 = (-t151 * t170 + t153 * t166) / t143 ^ 2;
	t206 = t129 * t145;
	t204 = t134 * t144;
	t203 = t137 * t158;
	t200 = t148 * t150;
	t197 = t155 * t156;
	t186 = 0.2e1 * t212;
	t185 = 0.2e1 * t211;
	t184 = t155 * t213;
	t183 = t123 * t197;
	t182 = t155 * t207;
	t181 = t140 * t200;
	t177 = 0.1e1 + t199;
	t176 = t207 * t219;
	t175 = 0.2e1 * t134 * t210;
	t172 = t158 * t181;
	t171 = t214 * t155 * t137;
	t169 = t177 * t156;
	t168 = t158 * t174;
	t167 = t130 * t204 + t206;
	t133 = -t157 * t194 - t196;
	t132 = -t156 * t144 + t180;
	t127 = t177 * t158 * t140;
	t118 = 0.1e1 / t121;
	t116 = 0.1e1 / t120;
	t112 = (-t138 * t172 + t171) * t156;
	t111 = t157 * t203 - t202 + (-t137 * t157 + t138 * t192) * t127;
	t110 = -t177 * t176 + (-qJD(1) * t169 + t166 * t219) * t140;
	t106 = -0.2e1 * t212 + 0.2e1 * (-t113 * t130 * t118 + (-t118 * t210 - t130 * t212) * t134) * t134;
	t1 = [t150 * t182 * t215 + (-qJD(3) * t169 - t190 * t198) * t140, 0, t110, 0, 0, 0; (-0.2e1 * t122 * t184 + (t122 * t188 + (-qJD(1) * t112 - t109) * t208) * t116) * t158 + (0.2e1 * t123 * t184 * t112 + (-((t115 * t172 + t214 * t188 - 0.2e1 * t182) * t137 + (t176 * t200 - t209 + (t209 + (-0.2e1 * t155 - t218) * t187) * t140) * t138) * t183 + (-t123 * t188 + t155 * t185) * t112 + (-t122 + ((-t149 + t153) * t138 * t181 - t158 * t171) * t123) * t155 * qJD(1)) * t116) * t156, 0, (t111 * t208 - t122 * t157) * t213 * t215 + ((t122 * t190 + (-qJD(3) * t111 - t109) * t156 * t123) * t157 + (-t122 * t189 - (t110 * t138 * t158 + t217 * t137 + (qJD(3) * t137 - t115 * t203 - t138 * t191) * t127) * t183 + (-t123 * t190 + t156 * t185) * t111 - ((-t110 - t191) * t137 + ((t127 * t158 - 0.1e1) * qJD(3) + (-t127 + t158) * t115) * t138) * t123 * t195) * t155) * t116, 0, 0, 0; (-t129 * t132 + t133 * t205) * t186 + (t133 * t175 - t129 * t144 * t168 - t216 * t206 + (t134 * t145 * t168 + t133 * t113 - t132 * t114 - t216 * t204) * t130) * t118, 0, t167 * t186 * t197 + (-t167 * t156 * t188 + (-t167 * t190 + ((t129 * t146 + t175) * t144 + (t113 * t144 + (-t134 * t146 + t114) * t145) * t130) * t156) * t155) * t118, 0, t106, t106;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end