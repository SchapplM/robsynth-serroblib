% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR3
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
%   Wie in S6RPRPRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:49
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:42
	% EndTime: 2019-10-10 00:49:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:42
	% EndTime: 2019-10-10 00:49:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:42
	% EndTime: 2019-10-10 00:49:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:42
	% EndTime: 2019-10-10 00:49:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:43
	% EndTime: 2019-10-10 00:49:44
	% DurationCPUTime: 0.92s
	% Computational Cost: add. (1727->84), mult. (2191->201), div. (456->12), fcn. (2616->9), ass. (0->86)
	t108 = sin(qJ(3));
	t102 = t108 ^ 2;
	t101 = t108 * t102;
	t109 = cos(qJ(3));
	t103 = 0.1e1 / t109;
	t104 = 0.1e1 / t109 ^ 2;
	t117 = qJD(3) * (t101 * t104 + t108) * t103;
	t100 = qJ(1) + pkin(10);
	t99 = cos(t100);
	t142 = qJD(1) * t99;
	t98 = sin(t100);
	t130 = t98 * t142;
	t140 = t102 * t104;
	t96 = t98 ^ 2;
	t95 = t96 * t140 + 0.1e1;
	t151 = (t96 * t117 + t130 * t140) / t95 ^ 2;
	t161 = -0.2e1 * t151;
	t106 = sin(pkin(11));
	t107 = cos(pkin(11));
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
	% StartTime: 2019-10-10 00:49:43
	% EndTime: 2019-10-10 00:49:44
	% DurationCPUTime: 1.10s
	% Computational Cost: add. (2324->95), mult. (2519->212), div. (480->12), fcn. (2968->9), ass. (0->95)
	t132 = qJ(1) + pkin(10);
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
	t131 = pkin(11) + qJ(5);
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
	% StartTime: 2019-10-10 00:49:43
	% EndTime: 2019-10-10 00:49:44
	% DurationCPUTime: 1.10s
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
	t150 = pkin(11) + qJ(5) + qJ(6);
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
	t151 = qJD(5) + qJD(6);
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
	t1 = [t149 * t155 * t177 + (-t148 * t155 * t194 + qJD(3) * t172) * t140, 0, t109, 0, 0, 0; (t121 * t179 + (-t121 * t191 + (qJD(1) * t111 + t108) * t212) * t118) * t148 + (t122 * t179 * t111 + (-((t114 * t174 + t218 * t191 + t177) * t136 + (t178 * t198 - t213 + (t213 + (t219 - t220) * t193) * t140) * t137) * t186 + (-t122 * t191 + t158 * t188) * t111 + (-t121 + ((-t146 + t147) * t137 * t184 + t218 * t185) * t122) * t194) * t118) * t149, 0, (t110 * t212 - t121 * t159) * t149 * t190 + ((-t121 * t196 + (-qJD(3) * t110 - t108) * t149 * t122) * t159 + (-t149 * qJD(3) * t121 - (-t109 * t137 * t148 + t136 * t193 + t210 * t214 - t214 + (-qJD(3) * t136 - t137 * t195) * t126) * t186 + (t122 * t196 + t149 * t188) * t110 - ((t109 - t195) * t136 + ((0.1e1 - t210) * qJD(3) + (t126 - t148) * t114) * t137) * t122 * t199) * t158) * t118, 0, 0, 0; (t128 * t170 + t132 * t208) * t189 + (0.2e1 * t132 * t187 - t176 * t128 * t203 + t167 * t209 + (-t176 * t133 * t204 + t132 * t112 + t113 * t170 - t167 * t207) * t129) * t115, 0, -t149 * t168 * t189 + (-t168 * t196 + (t171 * t191 + ((-t128 * t151 - 0.2e1 * t187) * t145 + (-t112 * t145 + (-t133 * t151 + t113) * t144) * t129) * t158) * t149) * t115, 0, t105, t105;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end