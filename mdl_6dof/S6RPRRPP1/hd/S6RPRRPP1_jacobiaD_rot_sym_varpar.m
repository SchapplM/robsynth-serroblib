% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPP1
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
%   Wie in S6RPRRPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:09
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:42
	% EndTime: 2019-10-10 01:09:43
	% DurationCPUTime: 1.02s
	% Computational Cost: add. (1976->91), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->94)
	t123 = qJ(1) + pkin(9);
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
	% StartTime: 2019-10-10 01:09:42
	% EndTime: 2019-10-10 01:09:43
	% DurationCPUTime: 1.03s
	% Computational Cost: add. (2324->95), mult. (2519->213), div. (480->12), fcn. (2968->9), ass. (0->94)
	t135 = qJ(1) + pkin(9);
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
	t134 = qJ(4) + pkin(10);
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
	% StartTime: 2019-10-10 01:09:42
	% EndTime: 2019-10-10 01:09:43
	% DurationCPUTime: 1.52s
	% Computational Cost: add. (7031->125), mult. (6168->274), div. (1114->15), fcn. (7752->9), ass. (0->114)
	t164 = qJ(4) + pkin(10);
	t161 = sin(t164);
	t162 = cos(t164);
	t216 = qJ(1) + pkin(9);
	t198 = sin(t216);
	t163 = cos(t216);
	t170 = cos(qJ(3));
	t224 = t163 * t170;
	t146 = t198 * t161 + t162 * t224;
	t140 = 0.1e1 / t146 ^ 2;
	t160 = t163 ^ 2;
	t169 = sin(qJ(3));
	t165 = t169 ^ 2;
	t227 = t160 * t165;
	t207 = t140 * t227;
	t135 = 0.1e1 + t207;
	t189 = qJD(1) * t198;
	t221 = qJD(3) * t163;
	t203 = t169 * t221;
	t178 = t170 * t189 + t203;
	t188 = t198 * qJD(4);
	t226 = t163 * t161;
	t125 = (-qJD(4) * t170 + qJD(1)) * t226 + (t188 - t178) * t162;
	t139 = 0.1e1 / t146;
	t237 = t125 * t139 * t140;
	t193 = t227 * t237;
	t220 = qJD(3) * t170;
	t201 = t169 * t220;
	t244 = (-t193 + (-t163 * t165 * t189 + t160 * t201) * t140) / t135 ^ 2;
	t225 = t163 * t169;
	t191 = t198 * t170;
	t142 = t161 * t191 + t163 * t162;
	t185 = t161 * t188;
	t218 = qJD(4) * t163;
	t199 = t162 * t218;
	t124 = t142 * qJD(1) + t161 * t203 - t170 * t199 - t185;
	t145 = t161 * t224 - t198 * t162;
	t157 = 0.1e1 / t161;
	t158 = 0.1e1 / t161 ^ 2;
	t166 = 0.1e1 / t169;
	t167 = 0.1e1 / t169 ^ 2;
	t202 = t167 * t220;
	t219 = qJD(4) * t162;
	t230 = t157 * t166;
	t243 = (t158 * t166 * t219 + t157 * t202) * t145 + t124 * t230;
	t223 = t169 * t161;
	t134 = atan2(-t142, t223);
	t129 = cos(t134);
	t128 = sin(t134);
	t236 = t128 * t142;
	t123 = t129 * t223 - t236;
	t120 = 0.1e1 / t123;
	t121 = 0.1e1 / t123 ^ 2;
	t242 = 0.2e1 * t145;
	t137 = t142 ^ 2;
	t228 = t158 * t167;
	t136 = t137 * t228 + 0.1e1;
	t132 = 0.1e1 / t136;
	t217 = qJD(4) * t169;
	t183 = t161 * t220 + t162 * t217;
	t205 = t142 * t228;
	t192 = t198 * t169;
	t186 = qJD(3) * t192;
	t187 = t162 * t189;
	t222 = qJD(1) * t163;
	t126 = t162 * t188 * t170 - t187 + (t222 * t170 - t186 - t218) * t161;
	t208 = t126 * t230;
	t112 = (t183 * t205 - t208) * t132;
	t179 = -t112 * t142 + t183;
	t108 = (-t112 * t223 - t126) * t128 + t179 * t129;
	t122 = t120 * t121;
	t241 = t108 * t122;
	t159 = t157 * t158;
	t168 = t166 / t165;
	t200 = t167 * t219;
	t240 = (t126 * t205 + (-t158 * t168 * t220 - t159 * t200) * t137) / t136 ^ 2;
	t239 = t121 * t145;
	t238 = t124 * t121;
	t235 = t128 * t145;
	t234 = t128 * t169;
	t233 = t129 * t142;
	t232 = t129 * t145;
	t231 = t129 * t170;
	t229 = t158 * t162;
	t138 = t145 ^ 2;
	t118 = t121 * t138 + 0.1e1;
	t215 = 0.2e1 * (-t138 * t241 - t145 * t238) / t118 ^ 2;
	t214 = 0.2e1 * t244;
	t213 = -0.2e1 * t240;
	t212 = t122 * t242;
	t211 = t166 * t240;
	t210 = t121 * t235;
	t206 = t142 * t230;
	t204 = t157 * t167 * t170;
	t197 = t120 * t215;
	t196 = t121 * t215;
	t195 = t225 * t242;
	t194 = t157 * t211;
	t182 = t142 * t204 + t198;
	t119 = t182 * t132;
	t190 = t198 - t119;
	t144 = t162 * t191 - t226;
	t184 = t142 * t229 - t144 * t157;
	t181 = t140 * t144 * t163 - t198 * t139;
	t130 = 0.1e1 / t135;
	t127 = t146 * qJD(1) - t162 * t186 - t170 * t185 - t199;
	t116 = 0.1e1 / t118;
	t115 = t184 * t166 * t132;
	t111 = (-t128 + (t129 * t206 + t128) * t132) * t145;
	t110 = -t119 * t233 + (t190 * t234 + t231) * t161;
	t109 = t129 * t162 * t169 - t128 * t144 + (-t128 * t223 - t233) * t115;
	t107 = t182 * t213 + (t126 * t204 + t222 + (-t158 * t170 * t200 + (-0.2e1 * t168 * t170 ^ 2 - t166) * t157 * qJD(3)) * t142) * t132;
	t105 = -0.2e1 * t184 * t211 + (-t184 * t202 + (t126 * t229 - t127 * t157 + (t144 * t229 + (-0.2e1 * t159 * t162 ^ 2 - t157) * t142) * qJD(4)) * t166) * t132;
	t1 = [t243 * t132 + t194 * t242, 0, t107, t105, 0, 0; t142 * t197 + (-t126 * t120 + (t108 * t142 + t111 * t124) * t121) * t116 + (t111 * t196 + (0.2e1 * t111 * t241 + (t124 * t132 - t124 - (-t112 * t132 * t206 + t213) * t145) * t121 * t128 + (-(-0.2e1 * t142 * t194 - t112) * t239 + (-(t112 + t208) * t145 + t243 * t142) * t121 * t132) * t129) * t116) * t145, 0, t110 * t145 * t196 + (-(-t107 * t233 + (t112 * t236 - t126 * t129) * t119) * t239 + (t108 * t212 + t238) * t110 + (-t120 * t225 - (-t119 * t234 + t128 * t192 + t231) * t239) * t219) * t116 + (t197 * t225 + ((-t120 * t221 - (t190 * qJD(3) - t112) * t210) * t170 + (t120 * t189 + (t163 * t108 - (-t107 + t222) * t235 - (t190 * t112 - qJD(3)) * t232) * t121) * t169) * t116) * t161, (t109 * t239 - t120 * t146) * t215 + (t109 * t238 + t125 * t120 + (t109 * t212 - t121 * t146) * t108 - (t162 * t220 - t161 * t217 - t105 * t142 - t115 * t126 + (-t115 * t223 - t144) * t112) * t121 * t232 - (-t127 + (-t105 * t161 - t112 * t162) * t169 - t179 * t115) * t210) * t116, 0, 0; t181 * t169 * t214 + (-t181 * t220 + ((qJD(1) * t139 + 0.2e1 * t144 * t237) * t163 + (-t198 * t125 - t127 * t163 + t144 * t189) * t140) * t169) * t130, 0, (t139 * t224 + t162 * t207) * t214 + (0.2e1 * t162 * t193 + t178 * t139 + ((t125 * t170 + 0.2e1 * t165 * t187) * t163 + (qJD(4) * t161 * t165 - 0.2e1 * t162 * t201) * t160) * t140) * t130, t140 * t195 * t244 + (t195 * t237 + (t124 * t225 + (-t163 * t220 + t169 * t189) * t145) * t140) * t130, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end