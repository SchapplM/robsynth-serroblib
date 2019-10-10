% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRP3
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
%   Wie in S6RPRPRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:32
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:32:30
	% EndTime: 2019-10-10 00:32:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:32:30
	% EndTime: 2019-10-10 00:32:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:32:30
	% EndTime: 2019-10-10 00:32:30
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:32:30
	% EndTime: 2019-10-10 00:32:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:32:31
	% EndTime: 2019-10-10 00:32:31
	% DurationCPUTime: 0.92s
	% Computational Cost: add. (1727->84), mult. (2191->201), div. (456->12), fcn. (2616->9), ass. (0->86)
	t108 = sin(qJ(3));
	t102 = t108 ^ 2;
	t101 = t108 * t102;
	t109 = cos(qJ(3));
	t103 = 0.1e1 / t109;
	t104 = 0.1e1 / t109 ^ 2;
	t117 = qJD(3) * (t101 * t104 + t108) * t103;
	t100 = qJ(1) + pkin(9);
	t99 = cos(t100);
	t142 = qJD(1) * t99;
	t98 = sin(t100);
	t130 = t98 * t142;
	t140 = t102 * t104;
	t96 = t98 ^ 2;
	t95 = t96 * t140 + 0.1e1;
	t151 = (t96 * t117 + t130 * t140) / t95 ^ 2;
	t161 = -0.2e1 * t151;
	t106 = sin(pkin(10));
	t107 = cos(pkin(10));
	t138 = t107 * t109;
	t89 = t98 * t106 + t99 * t138;
	t84 = 0.1e1 / t89 ^ 2;
	t139 = t106 * t109;
	t88 = -t98 * t107 + t99 * t139;
	t148 = t84 * t88;
	t136 = qJD(3) * t108;
	t124 = t107 * t136;
	t87 = t99 * t106 - t98 * t138;
	t81 = t87 * qJD(1) - t99 * t124;
	t159 = t81 * t84;
	t83 = 0.1e1 / t89;
	t149 = t83 * t159;
	t82 = t88 ^ 2;
	t77 = t82 * t84 + 0.1e1;
	t125 = t106 * t136;
	t86 = -t99 * t107 - t98 * t139;
	t80 = t86 * qJD(1) - t99 * t125;
	t160 = (t80 * t148 - t82 * t149) / t77 ^ 2;
	t122 = 0.1e1 + t140;
	t158 = t122 * t98;
	t144 = t98 * t108;
	t92 = atan2(-t144, -t109);
	t90 = sin(t92);
	t131 = t90 * t144;
	t91 = cos(t92);
	t76 = -t109 * t91 - t131;
	t71 = 0.1e1 / t76;
	t157 = -0.2e1 * t99;
	t72 = 0.1e1 / t76 ^ 2;
	t93 = 0.1e1 / t95;
	t156 = t93 - 0.1e1;
	t135 = qJD(3) * t109;
	t126 = t72 * t135;
	t141 = qJD(3) * t98;
	t145 = t109 * t90;
	t127 = t104 * t141;
	t137 = qJD(1) * t108;
	t67 = (-(-t98 * t135 - t99 * t137) * t103 + t102 * t127) * t93;
	t62 = (t67 - t141) * t145 + (-t90 * t142 + (-t67 * t98 + qJD(3)) * t91) * t108;
	t154 = t62 * t71 * t72;
	t97 = t99 ^ 2;
	t70 = t102 * t72 * t97 + 0.1e1;
	t155 = (t108 * t97 * t126 + (-t72 * t130 - t97 * t154) * t102) / t70 ^ 2;
	t153 = t67 * t90;
	t152 = t72 * t99;
	t79 = t93 * t158;
	t150 = t79 * t98;
	t147 = t79 - t98;
	t146 = t108 * t99;
	t143 = qJD(1) * t98;
	t134 = 0.2e1 * t154;
	t133 = t71 * t155;
	t132 = t88 * t149;
	t129 = t102 * t103 * t93;
	t128 = 0.1e1 - t150;
	t123 = 0.2e1 * t72 * t155;
	t121 = t103 * t161;
	t120 = t98 * t129;
	t119 = t122 * t99;
	t118 = -t106 * t83 + t107 * t148;
	t74 = 0.1e1 / t77;
	t68 = 0.1e1 / t70;
	t66 = (t156 * t90 * t108 - t91 * t120) * t99;
	t64 = t128 * t91 * t108 + t147 * t145;
	t63 = t158 * t161 + (qJD(1) * t119 + 0.2e1 * t117 * t98) * t93;
	t1 = [t121 * t146 + (-t103 * t98 * t137 + qJD(3) * t119) * t93, 0, t63, 0, 0, 0; (-t71 * t68 * t135 + (0.2e1 * t133 + (qJD(1) * t66 + t62) * t72 * t68) * t108) * t98 + (t66 * t123 * t108 + (-t66 * t126 + (t66 * t134 + ((0.2e1 * t108 * t151 - t67 * t120 - t156 * t135) * t90 + (t102 * t98 * t121 + t108 * t67 + (t101 * t127 - (t67 - 0.2e1 * t141) * t108) * t93) * t91) * t152) * t108 + (-t71 + (t156 * t131 - (t96 - t97) * t91 * t129) * t72) * t137) * t68) * t99, 0, (t133 * t157 + (-t71 * t143 + (-qJD(3) * t64 - t62) * t152) * t68) * t109 + (t64 * t99 * t123 + (-t99 * qJD(3) * t71 - (-t63 * t91 * t98 + t90 * t141 + t150 * t153 - t153 + (-qJD(3) * t90 - t142 * t91) * t79) * t72 * t146 + (t99 * t134 + t72 * t143) * t64 - ((t63 - t142) * t90 + (t128 * qJD(3) + t147 * t67) * t91) * t109 * t152) * t68) * t108, 0, 0, 0; 0.2e1 * (t87 * t148 - t83 * t86) * t160 + ((-t88 * qJD(1) + t98 * t125) * t83 + 0.2e1 * t87 * t132 + (-t86 * t81 - (-t89 * qJD(1) + t98 * t124) * t88 - t87 * t80) * t84) * t74, 0, t118 * t99 * t74 * t135 + (t118 * t157 * t160 + ((t83 * t143 + t99 * t159) * t106 + (t132 * t157 + (-t88 * t143 + t80 * t99) * t84) * t107) * t74) * t108, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:32:31
	% EndTime: 2019-10-10 00:32:32
	% DurationCPUTime: 1.03s
	% Computational Cost: add. (2324->95), mult. (2519->212), div. (480->12), fcn. (2968->9), ass. (0->95)
	t132 = qJ(1) + pkin(9);
	t128 = sin(t132);
	t125 = t128 ^ 2;
	t138 = sin(qJ(3));
	t134 = t138 ^ 2;
	t139 = cos(qJ(3));
	t136 = 0.1e1 / t139 ^ 2;
	t179 = t134 * t136;
	t122 = t125 * t179 + 0.1e1;
	t133 = t138 * t134;
	t135 = 0.1e1 / t139;
	t147 = qJD(3) * (t133 * t136 + t138) * t135;
	t130 = cos(t132);
	t177 = qJD(1) * t130;
	t185 = t128 * t134;
	t152 = t177 * t185;
	t193 = (t125 * t147 + t136 * t152) / t122 ^ 2;
	t202 = -0.2e1 * t193;
	t131 = pkin(10) + qJ(5);
	t127 = sin(t131);
	t129 = cos(t131);
	t180 = t130 * t139;
	t115 = t128 * t127 + t129 * t180;
	t110 = 0.1e1 / t115 ^ 2;
	t186 = t128 * t129;
	t114 = t127 * t180 - t186;
	t190 = t114 * t129;
	t109 = 0.1e1 / t115;
	t192 = t109 * t127;
	t149 = t110 * t190 - t192;
	t108 = t114 ^ 2;
	t101 = t108 * t110 + 0.1e1;
	t99 = 0.1e1 / t101;
	t201 = t149 * t99;
	t159 = 0.1e1 + t179;
	t200 = t128 * t159;
	t184 = t128 * t138;
	t119 = atan2(-t184, -t139);
	t117 = cos(t119);
	t116 = sin(t119);
	t166 = t116 * t184;
	t105 = -t117 * t139 - t166;
	t102 = 0.1e1 / t105;
	t103 = 0.1e1 / t105 ^ 2;
	t120 = 0.1e1 / t122;
	t199 = t120 - 0.1e1;
	t126 = t130 ^ 2;
	t173 = qJD(3) * t139;
	t188 = t126 * t134;
	t175 = qJD(3) * t128;
	t161 = t136 * t175;
	t176 = qJD(1) * t138;
	t162 = t130 * t176;
	t95 = (-(-t128 * t173 - t162) * t135 + t134 * t161) * t120;
	t156 = t95 - t175;
	t157 = -t128 * t95 + qJD(3);
	t189 = t117 * t138;
	t89 = t157 * t189 + (t156 * t139 - t162) * t116;
	t195 = t102 * t103 * t89;
	t98 = t103 * t188 + 0.1e1;
	t198 = (-t188 * t195 + (t126 * t138 * t173 - t152) * t103) / t98 ^ 2;
	t191 = t110 * t114;
	t154 = -qJD(1) * t139 + qJD(5);
	t155 = qJD(5) * t139 - qJD(1);
	t174 = qJD(3) * t138;
	t160 = t130 * t174;
	t182 = t130 * t127;
	t94 = -t155 * t182 + (t154 * t128 - t160) * t129;
	t194 = t109 * t110 * t94;
	t183 = t128 * t139;
	t148 = t127 * t183 + t130 * t129;
	t93 = t148 * qJD(1) - t115 * qJD(5) + t127 * t160;
	t197 = 0.1e1 / t101 ^ 2 * (-t108 * t194 - t93 * t191);
	t96 = 0.1e1 / t98;
	t196 = t103 * t96;
	t181 = t130 * t138;
	t178 = qJD(1) * t128;
	t172 = 0.2e1 * t198;
	t171 = -0.2e1 * t197;
	t170 = 0.2e1 * t195;
	t169 = t102 * t198;
	t168 = t114 * t194;
	t167 = t96 * t173;
	t165 = t120 * t134 * t135;
	t163 = t128 * t176;
	t158 = t135 * t202;
	t153 = t128 * t165;
	t151 = t159 * t130;
	t150 = t154 * t130;
	t113 = -t129 * t183 + t182;
	t107 = t120 * t200;
	t92 = (t199 * t138 * t116 - t117 * t153) * t130;
	t91 = -t116 * t183 + t189 + (t116 * t139 - t117 * t184) * t107;
	t90 = t200 * t202 + (qJD(1) * t151 + 0.2e1 * t128 * t147) * t120;
	t1 = [t158 * t181 + (qJD(3) * t151 - t135 * t163) * t120, 0, t90, 0, 0, 0; (-t102 * t167 + (0.2e1 * t169 + (qJD(1) * t92 + t89) * t196) * t138) * t128 + (t92 * t138 * t96 * t170 + (-t92 * t167 + (t92 * t172 + ((0.2e1 * t138 * t193 - t95 * t153 - t199 * t173) * t116 + (t158 * t185 + t138 * t95 + (t133 * t161 - (t95 - 0.2e1 * t175) * t138) * t120) * t117) * t96 * t130) * t138) * t103 + (-t102 + ((-t125 + t126) * t117 * t165 + t199 * t166) * t103) * t96 * t176) * t130, 0, (-t102 * t96 * t178 + (-0.2e1 * t169 + (-qJD(3) * t91 - t89) * t196) * t130) * t139 + (t91 * t130 * t103 * t172 + ((-qJD(3) * t102 + t91 * t170) * t130 + (t91 * t178 + (-(-t107 * t177 - t128 * t90) * t117 - ((t107 * t128 - 0.1e1) * t95 + (-t107 + t128) * qJD(3)) * t116) * t181) * t103) * t96 - ((t90 - t177) * t116 + (t156 * t107 + t157) * t117) * t180 * t196) * t138, 0, 0, 0; 0.2e1 * (t109 * t148 + t113 * t191) * t197 + (0.2e1 * t113 * t168 - t155 * t109 * t186 + (t128 * t174 + t150) * t192 + (t148 * t94 + t113 * t93 - t150 * t190 - (t155 * t127 + t129 * t174) * t114 * t128) * t110) * t99, 0, -t163 * t201 + (t173 * t201 + (t149 * t171 + ((-qJD(5) * t109 - 0.2e1 * t168) * t129 + (-t129 * t93 + (-qJD(5) * t114 + t94) * t127) * t110) * t99) * t138) * t130, 0, t171 + 0.2e1 * (-t110 * t93 * t99 + (-t110 * t197 - t99 * t194) * t114) * t114, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:32:31
	% EndTime: 2019-10-10 00:32:32
	% DurationCPUTime: 1.53s
	% Computational Cost: add. (7031->125), mult. (6168->274), div. (1114->15), fcn. (7752->9), ass. (0->114)
	t161 = pkin(10) + qJ(5);
	t158 = sin(t161);
	t159 = cos(t161);
	t213 = qJ(1) + pkin(9);
	t195 = sin(t213);
	t160 = cos(t213);
	t167 = cos(qJ(3));
	t221 = t160 * t167;
	t143 = t195 * t158 + t159 * t221;
	t137 = 0.1e1 / t143 ^ 2;
	t157 = t160 ^ 2;
	t166 = sin(qJ(3));
	t162 = t166 ^ 2;
	t224 = t157 * t162;
	t204 = t137 * t224;
	t132 = 0.1e1 + t204;
	t186 = qJD(1) * t195;
	t218 = qJD(3) * t160;
	t200 = t166 * t218;
	t175 = t167 * t186 + t200;
	t185 = t195 * qJD(5);
	t223 = t160 * t158;
	t122 = (-qJD(5) * t167 + qJD(1)) * t223 + (t185 - t175) * t159;
	t136 = 0.1e1 / t143;
	t234 = t122 * t136 * t137;
	t190 = t224 * t234;
	t217 = qJD(3) * t167;
	t198 = t166 * t217;
	t241 = (-t190 + (-t160 * t162 * t186 + t157 * t198) * t137) / t132 ^ 2;
	t222 = t160 * t166;
	t188 = t195 * t167;
	t139 = t158 * t188 + t160 * t159;
	t182 = t158 * t185;
	t215 = qJD(5) * t160;
	t196 = t159 * t215;
	t121 = t139 * qJD(1) + t158 * t200 - t167 * t196 - t182;
	t142 = t158 * t221 - t195 * t159;
	t154 = 0.1e1 / t158;
	t155 = 0.1e1 / t158 ^ 2;
	t163 = 0.1e1 / t166;
	t164 = 0.1e1 / t166 ^ 2;
	t199 = t164 * t217;
	t216 = qJD(5) * t159;
	t227 = t154 * t163;
	t240 = (t155 * t163 * t216 + t154 * t199) * t142 + t121 * t227;
	t220 = t166 * t158;
	t131 = atan2(-t139, t220);
	t126 = cos(t131);
	t125 = sin(t131);
	t233 = t125 * t139;
	t120 = t126 * t220 - t233;
	t117 = 0.1e1 / t120;
	t118 = 0.1e1 / t120 ^ 2;
	t239 = 0.2e1 * t142;
	t134 = t139 ^ 2;
	t225 = t155 * t164;
	t133 = t134 * t225 + 0.1e1;
	t129 = 0.1e1 / t133;
	t214 = qJD(5) * t166;
	t180 = t158 * t217 + t159 * t214;
	t202 = t139 * t225;
	t189 = t195 * t166;
	t183 = qJD(3) * t189;
	t184 = t159 * t186;
	t219 = qJD(1) * t160;
	t123 = t159 * t185 * t167 - t184 + (t219 * t167 - t183 - t215) * t158;
	t205 = t123 * t227;
	t109 = (t180 * t202 - t205) * t129;
	t176 = -t109 * t139 + t180;
	t105 = (-t109 * t220 - t123) * t125 + t176 * t126;
	t119 = t117 * t118;
	t238 = t105 * t119;
	t156 = t154 * t155;
	t165 = t163 / t162;
	t197 = t164 * t216;
	t237 = (t123 * t202 + (-t155 * t165 * t217 - t156 * t197) * t134) / t133 ^ 2;
	t236 = t118 * t142;
	t235 = t121 * t118;
	t232 = t125 * t142;
	t231 = t125 * t166;
	t230 = t126 * t139;
	t229 = t126 * t142;
	t228 = t126 * t167;
	t226 = t155 * t159;
	t135 = t142 ^ 2;
	t115 = t118 * t135 + 0.1e1;
	t212 = 0.2e1 * (-t135 * t238 - t142 * t235) / t115 ^ 2;
	t211 = 0.2e1 * t241;
	t210 = -0.2e1 * t237;
	t209 = t119 * t239;
	t208 = t163 * t237;
	t207 = t118 * t232;
	t203 = t139 * t227;
	t201 = t154 * t164 * t167;
	t194 = t117 * t212;
	t193 = t118 * t212;
	t192 = t222 * t239;
	t191 = t154 * t208;
	t179 = t139 * t201 + t195;
	t116 = t179 * t129;
	t187 = t195 - t116;
	t141 = t159 * t188 - t223;
	t181 = t139 * t226 - t141 * t154;
	t178 = t137 * t141 * t160 - t195 * t136;
	t127 = 0.1e1 / t132;
	t124 = t143 * qJD(1) - t159 * t183 - t167 * t182 - t196;
	t113 = 0.1e1 / t115;
	t112 = t181 * t163 * t129;
	t108 = (-t125 + (t126 * t203 + t125) * t129) * t142;
	t107 = -t116 * t230 + (t187 * t231 + t228) * t158;
	t106 = t126 * t159 * t166 - t125 * t141 + (-t125 * t220 - t230) * t112;
	t104 = t179 * t210 + (t123 * t201 + t219 + (-t155 * t167 * t197 + (-0.2e1 * t165 * t167 ^ 2 - t163) * t154 * qJD(3)) * t139) * t129;
	t102 = -0.2e1 * t181 * t208 + (-t181 * t199 + (t123 * t226 - t124 * t154 + (t141 * t226 + (-0.2e1 * t156 * t159 ^ 2 - t154) * t139) * qJD(5)) * t163) * t129;
	t1 = [t129 * t240 + t191 * t239, 0, t104, 0, t102, 0; t139 * t194 + (-t123 * t117 + (t105 * t139 + t108 * t121) * t118) * t113 + (t108 * t193 + (0.2e1 * t108 * t238 + (t121 * t129 - t121 - (-t109 * t129 * t203 + t210) * t142) * t118 * t125 + (-(-0.2e1 * t139 * t191 - t109) * t236 + (-(t109 + t205) * t142 + t240 * t139) * t118 * t129) * t126) * t113) * t142, 0, t107 * t142 * t193 + (-(-t104 * t230 + (t109 * t233 - t123 * t126) * t116) * t236 + (t105 * t209 + t235) * t107 + (-t117 * t222 - (-t116 * t231 + t125 * t189 + t228) * t236) * t216) * t113 + (t194 * t222 + ((-t117 * t218 - (t187 * qJD(3) - t109) * t207) * t167 + (t117 * t186 + (t160 * t105 - (-t104 + t219) * t232 - (t187 * t109 - qJD(3)) * t229) * t118) * t166) * t113) * t158, 0, (t106 * t236 - t117 * t143) * t212 + (t106 * t235 + t122 * t117 + (t106 * t209 - t118 * t143) * t105 - (t159 * t217 - t158 * t214 - t102 * t139 - t112 * t123 + (-t112 * t220 - t141) * t109) * t118 * t229 - (-t124 + (-t102 * t158 - t109 * t159) * t166 - t176 * t112) * t207) * t113, 0; t178 * t166 * t211 + (-t178 * t217 + ((qJD(1) * t136 + 0.2e1 * t141 * t234) * t160 + (-t195 * t122 - t124 * t160 + t141 * t186) * t137) * t166) * t127, 0, (t136 * t221 + t159 * t204) * t211 + (0.2e1 * t159 * t190 + t175 * t136 + ((t122 * t167 + 0.2e1 * t162 * t184) * t160 + (qJD(5) * t158 * t162 - 0.2e1 * t159 * t198) * t157) * t137) * t127, 0, t137 * t192 * t241 + (t192 * t234 + (t121 * t222 + (-t160 * t217 + t166 * t186) * t142) * t137) * t127, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end