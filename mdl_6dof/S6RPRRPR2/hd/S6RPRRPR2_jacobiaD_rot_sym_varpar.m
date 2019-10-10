% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR2
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
%   Wie in S6RPRRPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:15
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (1976->91), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->94)
	t123 = qJ(1) + pkin(10);
	t121 = sin(t123);
	t119 = t121 ^ 2;
	t130 = sin(qJ(3));
	t125 = t130 ^ 2;
	t132 = cos(qJ(3));
	t127 = 0.1e1 / t132 ^ 2;
	t175 = t125 * t127;
	t115 = t119 * t175 + 0.1e1;
	t124 = t130 * t125;
	t126 = 0.1e1 / t132;
	t140 = qJD(3) * (t124 * t127 + t130) * t126;
	t122 = cos(t123);
	t171 = qJD(1) * t122;
	t180 = t121 * t125;
	t145 = t171 * t180;
	t187 = 0.1e1 / t115 ^ 2 * (t119 * t140 + t127 * t145);
	t198 = -0.2e1 * t187;
	t129 = sin(qJ(4));
	t131 = cos(qJ(4));
	t173 = t131 * t132;
	t109 = t121 * t129 + t122 * t173;
	t103 = 0.1e1 / t109;
	t104 = 0.1e1 / t109 ^ 2;
	t174 = t129 * t132;
	t108 = -t121 * t131 + t122 * t174;
	t183 = t104 * t108;
	t142 = -t103 * t129 + t131 * t183;
	t102 = t108 ^ 2;
	t101 = t102 * t104 + 0.1e1;
	t98 = 0.1e1 / t101;
	t197 = t142 * t98;
	t169 = qJD(3) * t121;
	t113 = 0.1e1 / t115;
	t156 = t127 * t169;
	t170 = qJD(1) * t130;
	t157 = t122 * t170;
	t167 = qJD(3) * t132;
	t87 = (-(-t121 * t167 - t157) * t126 + t125 * t156) * t113;
	t149 = t87 - t169;
	t152 = 0.1e1 + t175;
	t196 = t121 * t152;
	t150 = -t121 * t87 + qJD(3);
	t148 = qJD(4) * t132 - qJD(1);
	t168 = qJD(3) * t130;
	t195 = t148 * t129 + t131 * t168;
	t178 = t121 * t130;
	t112 = atan2(-t178, -t132);
	t111 = cos(t112);
	t110 = sin(t112);
	t161 = t110 * t178;
	t96 = -t111 * t132 - t161;
	t93 = 0.1e1 / t96;
	t94 = 0.1e1 / t96 ^ 2;
	t194 = t113 - 0.1e1;
	t120 = t122 ^ 2;
	t162 = t94 * t167;
	t181 = t111 * t130;
	t82 = t150 * t181 + (t149 * t132 - t157) * t110;
	t192 = t82 * t93 * t94;
	t92 = t120 * t125 * t94 + 0.1e1;
	t193 = (-t94 * t145 + (-t125 * t192 + t130 * t162) * t120) / t92 ^ 2;
	t147 = -qJD(1) * t132 + qJD(4);
	t143 = t147 * t131;
	t89 = t121 * t143 - t195 * t122;
	t188 = t103 * t104 * t89;
	t141 = t121 * t174 + t122 * t131;
	t155 = t129 * t168;
	t88 = t141 * qJD(1) - t109 * qJD(4) + t122 * t155;
	t191 = (-t102 * t188 - t88 * t183) / t101 ^ 2;
	t90 = 0.1e1 / t92;
	t190 = t90 * t94;
	t189 = t93 * t90;
	t185 = t122 * t94;
	t182 = t110 * t132;
	t177 = t122 * t129;
	t176 = t122 * t130;
	t172 = qJD(1) * t121;
	t166 = 0.2e1 * t192;
	t165 = -0.2e1 * t191;
	t164 = t93 * t193;
	t163 = t108 * t188;
	t160 = t113 * t125 * t126;
	t158 = t121 * t170;
	t153 = 0.2e1 * t94 * t193;
	t151 = t126 * t198;
	t146 = t121 * t160;
	t144 = t152 * t122;
	t107 = -t121 * t173 + t177;
	t100 = t113 * t196;
	t86 = (t194 * t130 * t110 - t111 * t146) * t122;
	t85 = -t121 * t182 + t181 + (-t111 * t178 + t182) * t100;
	t83 = t196 * t198 + (qJD(1) * t144 + 0.2e1 * t121 * t140) * t113;
	t1 = [t151 * t176 + (qJD(3) * t144 - t126 * t158) * t113, 0, t83, 0, 0, 0; (-t167 * t189 + (0.2e1 * t164 + (qJD(1) * t86 + t82) * t190) * t130) * t121 + (t86 * t153 * t130 + (-t86 * t162 + (t86 * t166 + ((0.2e1 * t130 * t187 - t87 * t146 - t194 * t167) * t110 + (t151 * t180 + t130 * t87 + (t124 * t156 - (t87 - 0.2e1 * t169) * t130) * t113) * t111) * t185) * t130 + (-t93 + (-(t119 - t120) * t111 * t160 + t194 * t161) * t94) * t170) * t90) * t122, 0, (-t172 * t189 + (-0.2e1 * t164 + (-qJD(3) * t85 - t82) * t190) * t122) * t132 + (t85 * t122 * t153 + (-t122 * qJD(3) * t93 - ((-t100 * t171 - t121 * t83) * t111 + (-t150 * t100 - t149) * t110) * t94 * t176 + (t122 * t166 + t94 * t172) * t85 - ((t83 - t171) * t110 + (t149 * t100 + t150) * t111) * t132 * t185) * t90) * t130, 0, 0, 0; 0.2e1 * (t103 * t141 + t107 * t183) * t191 + (0.2e1 * t107 * t163 + (t141 * t89 + t107 * t88 + (-t195 * t121 - t122 * t143) * t108) * t104 + (t147 * t177 + (-t148 * t131 + t155) * t121) * t103) * t98, 0, -t158 * t197 + (t167 * t197 + (t142 * t165 + ((-qJD(4) * t103 - 0.2e1 * t163) * t131 + (-t131 * t88 + (-qJD(4) * t108 + t89) * t129) * t104) * t98) * t130) * t122, t165 + 0.2e1 * (-t104 * t88 * t98 + (-t104 * t191 - t98 * t188) * t108) * t108, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:15
	% DurationCPUTime: 1.04s
	% Computational Cost: add. (2324->95), mult. (2519->213), div. (480->12), fcn. (2968->9), ass. (0->94)
	t135 = qJ(1) + pkin(10);
	t131 = sin(t135);
	t128 = t131 ^ 2;
	t141 = sin(qJ(3));
	t137 = t141 ^ 2;
	t142 = cos(qJ(3));
	t139 = 0.1e1 / t142 ^ 2;
	t182 = t137 * t139;
	t125 = t128 * t182 + 0.1e1;
	t136 = t141 * t137;
	t138 = 0.1e1 / t142;
	t151 = qJD(3) * (t136 * t139 + t141) * t138;
	t133 = cos(t135);
	t180 = qJD(1) * t133;
	t188 = t131 * t137;
	t156 = t180 * t188;
	t196 = (t128 * t151 + t139 * t156) / t125 ^ 2;
	t204 = -0.2e1 * t196;
	t163 = 0.1e1 + t182;
	t203 = t131 * t163;
	t134 = qJ(4) + pkin(11);
	t130 = sin(t134);
	t132 = cos(t134);
	t183 = t133 * t142;
	t118 = t131 * t130 + t132 * t183;
	t187 = t131 * t141;
	t122 = atan2(-t187, -t142);
	t120 = cos(t122);
	t119 = sin(t122);
	t169 = t119 * t187;
	t108 = -t120 * t142 - t169;
	t105 = 0.1e1 / t108;
	t112 = 0.1e1 / t118;
	t106 = 0.1e1 / t108 ^ 2;
	t113 = 0.1e1 / t118 ^ 2;
	t123 = 0.1e1 / t125;
	t202 = t123 - 0.1e1;
	t129 = t133 ^ 2;
	t191 = t129 * t137;
	t101 = t106 * t191 + 0.1e1;
	t176 = qJD(3) * t142;
	t178 = qJD(3) * t131;
	t165 = t139 * t178;
	t179 = qJD(1) * t141;
	t166 = t133 * t179;
	t98 = (-(-t131 * t176 - t166) * t138 + t137 * t165) * t123;
	t160 = t98 - t178;
	t161 = -t131 * t98 + qJD(3);
	t192 = t120 * t141;
	t92 = t161 * t192 + (t160 * t142 - t166) * t119;
	t198 = t105 * t106 * t92;
	t201 = 0.1e1 / t101 ^ 2 * (-t191 * t198 + (t129 * t141 * t176 - t156) * t106);
	t189 = t131 * t132;
	t117 = t130 * t183 - t189;
	t111 = t117 ^ 2;
	t104 = t111 * t113 + 0.1e1;
	t194 = t113 * t117;
	t158 = -qJD(1) * t142 + qJD(4);
	t159 = qJD(4) * t142 - qJD(1);
	t177 = qJD(3) * t141;
	t164 = t133 * t177;
	t185 = t133 * t130;
	t97 = -t159 * t185 + (t158 * t131 - t164) * t132;
	t197 = t112 * t113 * t97;
	t186 = t131 * t142;
	t152 = t130 * t186 + t133 * t132;
	t96 = t152 * qJD(1) - qJD(4) * t118 + t130 * t164;
	t200 = 0.1e1 / t104 ^ 2 * (-t111 * t197 - t96 * t194);
	t99 = 0.1e1 / t101;
	t199 = t106 * t99;
	t195 = t112 * t130;
	t193 = t117 * t132;
	t184 = t133 * t141;
	t181 = qJD(1) * t131;
	t175 = 0.2e1 * t201;
	t174 = 0.2e1 * t200;
	t173 = 0.2e1 * t198;
	t172 = t105 * t201;
	t171 = t117 * t197;
	t170 = t99 * t176;
	t168 = t123 * t137 * t138;
	t162 = t138 * t204;
	t157 = t131 * t168;
	t155 = t163 * t133;
	t154 = t158 * t133;
	t153 = t113 * t193 - t195;
	t150 = t153 * t141;
	t116 = -t132 * t186 + t185;
	t110 = t123 * t203;
	t102 = 0.1e1 / t104;
	t95 = (t202 * t141 * t119 - t120 * t157) * t133;
	t94 = -t119 * t186 + t192 + (t119 * t142 - t120 * t187) * t110;
	t93 = t203 * t204 + (qJD(1) * t155 + 0.2e1 * t131 * t151) * t123;
	t1 = [t162 * t184 + (-t131 * t138 * t179 + qJD(3) * t155) * t123, 0, t93, 0, 0, 0; (-t105 * t170 + (0.2e1 * t172 + (qJD(1) * t95 + t92) * t199) * t141) * t131 + (t95 * t141 * t99 * t173 + (-t95 * t170 + (t95 * t175 + ((0.2e1 * t141 * t196 - t98 * t157 - t202 * t176) * t119 + (t162 * t188 + t141 * t98 + (t136 * t165 - (t98 - 0.2e1 * t178) * t141) * t123) * t120) * t99 * t133) * t141) * t106 + (-t105 + ((-t128 + t129) * t120 * t168 + t202 * t169) * t106) * t99 * t179) * t133, 0, (-t105 * t99 * t181 + (-0.2e1 * t172 + (-qJD(3) * t94 - t92) * t199) * t133) * t142 + (t94 * t133 * t106 * t175 + ((-qJD(3) * t105 + t94 * t173) * t133 + (t94 * t181 + (-(-t110 * t180 - t131 * t93) * t120 - ((t110 * t131 - 0.1e1) * t98 + (-t110 + t131) * qJD(3)) * t119) * t184) * t106) * t99 - ((t93 - t180) * t119 + (t160 * t110 + t161) * t120) * t183 * t199) * t141, 0, 0, 0; (t112 * t152 + t116 * t194) * t174 + (0.2e1 * t116 * t171 - t159 * t112 * t189 + (t131 * t177 + t154) * t195 + (t152 * t97 + t116 * t96 - t154 * t193 - (t159 * t130 + t132 * t177) * t117 * t131) * t113) * t102, 0, -t133 * t150 * t174 + (-t150 * t181 + (t153 * t176 + ((-qJD(4) * t112 - 0.2e1 * t171) * t132 + (-t132 * t96 + (-qJD(4) * t117 + t97) * t130) * t113) * t141) * t133) * t102, -0.2e1 * t200 + 0.2e1 * (-t102 * t113 * t96 + (-t102 * t197 - t113 * t200) * t117) * t117, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:15
	% DurationCPUTime: 1.11s
	% Computational Cost: add. (3074->97), mult. (2734->204), div. (498->12), fcn. (3199->9), ass. (0->96)
	t158 = sin(qJ(3));
	t154 = t158 ^ 2;
	t159 = cos(qJ(3));
	t156 = 0.1e1 / t159 ^ 2;
	t197 = t154 * t156;
	t152 = qJ(1) + pkin(10);
	t148 = sin(t152);
	t221 = 0.2e1 * t148;
	t220 = t158 * t197;
	t150 = qJ(4) + pkin(11) + qJ(6);
	t145 = cos(t150);
	t149 = cos(t152);
	t199 = t149 * t159;
	t144 = sin(t150);
	t204 = t148 * t144;
	t134 = t145 * t199 + t204;
	t202 = t148 * t158;
	t139 = atan2(-t202, -t159);
	t137 = cos(t139);
	t136 = sin(t139);
	t185 = t136 * t202;
	t124 = -t137 * t159 - t185;
	t121 = 0.1e1 / t124;
	t128 = 0.1e1 / t134;
	t155 = 0.1e1 / t159;
	t122 = 0.1e1 / t124 ^ 2;
	t129 = 0.1e1 / t134 ^ 2;
	t219 = -0.2e1 * t158;
	t146 = t148 ^ 2;
	t142 = t146 * t197 + 0.1e1;
	t140 = 0.1e1 / t142;
	t218 = t140 - 0.1e1;
	t151 = qJD(4) + qJD(6);
	t201 = t148 * t159;
	t170 = t144 * t201 + t149 * t145;
	t192 = qJD(3) * t158;
	t181 = t149 * t192;
	t112 = t170 * qJD(1) - t134 * t151 + t144 * t181;
	t203 = t148 * t145;
	t133 = t144 * t199 - t203;
	t127 = t133 ^ 2;
	t117 = t127 * t129 + 0.1e1;
	t208 = t129 * t133;
	t175 = -qJD(1) * t159 + t151;
	t176 = t151 * t159 - qJD(1);
	t200 = t149 * t144;
	t113 = -t176 * t200 + (t175 * t148 - t181) * t145;
	t215 = t113 * t128 * t129;
	t217 = (-t112 * t208 - t127 * t215) / t117 ^ 2;
	t194 = qJD(1) * t158;
	t182 = t149 * t194;
	t191 = qJD(3) * t159;
	t193 = qJD(3) * t148;
	t114 = (-(-t148 * t191 - t182) * t155 + t193 * t197) * t140;
	t206 = t137 * t158;
	t108 = (-t114 * t148 + qJD(3)) * t206 + (-t182 + (t114 - t193) * t159) * t136;
	t216 = t108 * t121 * t122;
	t214 = t114 * t136;
	t213 = t114 * t158;
	t212 = t122 * t158;
	t169 = qJD(3) * (t158 + t220) * t155;
	t195 = qJD(1) * t149;
	t173 = t148 * t154 * t195;
	t211 = (t146 * t169 + t156 * t173) / t142 ^ 2;
	t180 = 0.1e1 + t197;
	t126 = t180 * t148 * t140;
	t210 = t126 * t148;
	t209 = t128 * t144;
	t207 = t133 * t145;
	t147 = t149 ^ 2;
	t205 = t147 * t154;
	t198 = t154 * t155;
	t196 = qJD(1) * t148;
	t120 = t122 * t205 + 0.1e1;
	t190 = 0.2e1 * (-t205 * t216 + (t147 * t158 * t191 - t173) * t122) / t120 ^ 2;
	t189 = 0.2e1 * t217;
	t188 = 0.2e1 * t216;
	t187 = t133 * t215;
	t186 = t149 * t212;
	t184 = t140 * t198;
	t179 = t158 * t190;
	t178 = t211 * t221;
	t177 = t211 * t219;
	t174 = t148 * t184;
	t172 = t180 * t149;
	t171 = t129 * t207 - t209;
	t168 = t171 * t158;
	t167 = t148 * t192 + t175 * t149;
	t132 = -t145 * t201 + t200;
	t118 = 0.1e1 / t120;
	t115 = 0.1e1 / t117;
	t111 = (t218 * t158 * t136 - t137 * t174) * t149;
	t110 = -t136 * t201 + t206 + (t136 * t159 - t137 * t202) * t126;
	t109 = -t180 * t178 + (qJD(1) * t172 + t169 * t221) * t140;
	t105 = -0.2e1 * t217 + 0.2e1 * (-t112 * t115 * t129 + (-t115 * t215 - t129 * t217) * t133) * t133;
	t1 = [t149 * t155 * t177 + (-t148 * t155 * t194 + qJD(3) * t172) * t140, 0, t109, 0, 0, 0; (t121 * t179 + (-t121 * t191 + (qJD(1) * t111 + t108) * t212) * t118) * t148 + (t122 * t179 * t111 + (-((t114 * t174 + t218 * t191 + t177) * t136 + (t178 * t198 - t213 + (t213 + (t219 - t220) * t193) * t140) * t137) * t186 + (-t122 * t191 + t158 * t188) * t111 + (-t121 + ((-t146 + t147) * t137 * t184 + t218 * t185) * t122) * t194) * t118) * t149, 0, (t110 * t212 - t121 * t159) * t149 * t190 + ((-t121 * t196 + (-qJD(3) * t110 - t108) * t149 * t122) * t159 + (-t149 * qJD(3) * t121 - (-t109 * t137 * t148 + t136 * t193 + t210 * t214 - t214 + (-qJD(3) * t136 - t137 * t195) * t126) * t186 + (t122 * t196 + t149 * t188) * t110 - ((t109 - t195) * t136 + ((0.1e1 - t210) * qJD(3) + (t126 - t148) * t114) * t137) * t122 * t199) * t158) * t118, 0, 0, 0; (t128 * t170 + t132 * t208) * t189 + (0.2e1 * t132 * t187 - t176 * t128 * t203 + t167 * t209 + (-t176 * t133 * t204 + t132 * t112 + t113 * t170 - t167 * t207) * t129) * t115, 0, -t149 * t168 * t189 + (-t168 * t196 + (t171 * t191 + ((-t128 * t151 - 0.2e1 * t187) * t145 + (-t112 * t145 + (-t133 * t151 + t113) * t144) * t129) * t158) * t149) * t115, t105, 0, t105;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end