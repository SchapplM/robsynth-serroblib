% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRP9
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
%   Wie in S6RPRPRP9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:42
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP9_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP9_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:42:48
	% EndTime: 2019-10-10 00:42:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:42:48
	% EndTime: 2019-10-10 00:42:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:42:48
	% EndTime: 2019-10-10 00:42:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:42:48
	% EndTime: 2019-10-10 00:42:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:42:48
	% EndTime: 2019-10-10 00:42:49
	% DurationCPUTime: 0.96s
	% Computational Cost: add. (703->82), mult. (2191->190), div. (456->12), fcn. (2616->9), ass. (0->85)
	t102 = sin(qJ(3));
	t94 = 0.1e1 / t102 ^ 2;
	t104 = cos(qJ(3));
	t98 = t104 ^ 2;
	t147 = t94 * t98;
	t105 = cos(qJ(1));
	t126 = 0.1e1 + t147;
	t99 = t105 ^ 2;
	t92 = t99 * t147 + 0.1e1;
	t90 = 0.1e1 / t92;
	t115 = t126 * t90;
	t76 = t105 * t115;
	t161 = t105 * t76 - 0.1e1;
	t101 = cos(pkin(9));
	t103 = sin(qJ(1));
	t133 = qJD(3) * t104;
	t122 = t103 * t133;
	t100 = sin(pkin(9));
	t140 = t103 * t100;
	t141 = t102 * t105;
	t86 = t101 * t141 - t140;
	t78 = t86 * qJD(1) + t101 * t122;
	t139 = t103 * t101;
	t84 = t100 * t105 + t102 * t139;
	t81 = 0.1e1 / t84 ^ 2;
	t160 = t78 * t81;
	t159 = t104 * t147;
	t83 = -t101 * t105 + t102 * t140;
	t149 = t81 * t83;
	t79 = t83 ^ 2;
	t74 = t79 * t81 + 0.1e1;
	t72 = 0.1e1 / t74;
	t80 = 0.1e1 / t84;
	t158 = (-t100 * t80 + t101 * t149) * t72;
	t138 = t105 * t104;
	t89 = atan2(-t138, t102);
	t87 = sin(t89);
	t88 = cos(t89);
	t71 = t88 * t102 - t87 * t138;
	t68 = 0.1e1 / t71;
	t93 = 0.1e1 / t102;
	t69 = 0.1e1 / t71 ^ 2;
	t157 = t90 - 0.1e1;
	t135 = qJD(1) * t105;
	t116 = t103 * t98 * t135;
	t96 = t103 ^ 2;
	t146 = t96 * t98;
	t132 = qJD(3) * t105;
	t137 = qJD(1) * t103;
	t145 = t102 * t87;
	t136 = qJD(1) * t104;
	t64 = ((t102 * t132 + t103 * t136) * t93 + t132 * t147) * t90;
	t59 = (-t64 + t132) * t145 + (t87 * t137 + (-t105 * t64 + qJD(3)) * t88) * t104;
	t155 = t59 * t68 * t69;
	t67 = t69 * t146 + 0.1e1;
	t156 = (-t146 * t155 + (-t102 * t96 * t133 + t116) * t69) / t67 ^ 2;
	t150 = t80 * t160;
	t85 = t100 * t141 + t139;
	t77 = t85 * qJD(1) + t100 * t122;
	t154 = (t77 * t149 - t79 * t150) / t74 ^ 2;
	t65 = 0.1e1 / t67;
	t152 = t65 * t69;
	t113 = qJD(3) * (-t104 - t159) * t93;
	t151 = (t99 * t113 - t94 * t116) / t92 ^ 2;
	t148 = t93 * t98;
	t144 = t103 * t69;
	t142 = qJD(3) * t76;
	t134 = qJD(3) * t102;
	t131 = -0.2e1 * t155;
	t130 = 0.2e1 * t154;
	t129 = t68 * t156;
	t128 = t90 * t148;
	t127 = t104 * t151;
	t125 = t104 * t157;
	t124 = t65 * t134;
	t123 = t104 * t135;
	t121 = t104 * t132;
	t120 = -0.2e1 * t69 * t156;
	t119 = 0.2e1 * t83 * t150;
	t118 = t105 * t128;
	t117 = t87 * t125;
	t63 = (-t88 * t118 - t117) * t103;
	t61 = -t161 * t88 * t104 + (t105 - t76) * t145;
	t60 = -t115 * t137 + 0.2e1 * (t113 * t90 - t126 * t151) * t105;
	t1 = [t93 * t90 * t123 + (-qJD(3) * t115 - 0.2e1 * t93 * t127) * t103, 0, t60, 0, 0, 0; (t68 * t124 + (0.2e1 * t129 + (qJD(1) * t63 + t59) * t152) * t104) * t105 + (-t63 * t69 * t124 + (t63 * t120 + (t63 * t131 + ((t64 * t118 + t157 * t134 + 0.2e1 * t127) * t87 + (-t64 * t125 + (0.2e1 * t148 * t151 + (0.2e1 * t104 + t159) * t90 * qJD(3)) * t105) * t88) * t144) * t65) * t104 + (t68 + ((t96 - t99) * t88 * t128 - t105 * t117) * t69) * t65 * t136) * t103, 0, (t68 * t65 * t135 + (-0.2e1 * t129 + (-qJD(3) * t61 - t59) * t152) * t103) * t102 + (t61 * t103 * t120 + (t103 * qJD(3) * t68 + (t103 * t131 + t69 * t135) * t61 + (((-t105 * t60 + t137 * t76) * t88 + (t161 * t64 + t132 - t142) * t87) * t104 + ((-t60 - t137) * t87 + (-t64 * t76 - qJD(3) + (t64 + t142) * t105) * t88) * t102) * t144) * t65) * t104, 0, 0, 0; (t86 * t149 - t80 * t85) * t130 + ((-t83 * qJD(1) + t100 * t121) * t80 + t86 * t119 + (-t85 * t78 - (-t84 * qJD(1) + t101 * t121) * t83 - t86 * t77) * t81) * t72, 0, -t123 * t158 + (t134 * t158 + ((-0.2e1 * t80 * t154 - t72 * t160) * t100 + (t130 * t149 + (-t77 * t81 + t119) * t72) * t101) * t104) * t103, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:42:48
	% EndTime: 2019-10-10 00:42:49
	% DurationCPUTime: 0.94s
	% Computational Cost: add. (1161->91), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->95)
	t135 = cos(qJ(1));
	t198 = 0.2e1 * t135;
	t124 = pkin(9) + qJ(5);
	t122 = sin(t124);
	t132 = sin(qJ(3));
	t123 = cos(t124);
	t133 = sin(qJ(1));
	t179 = t133 * t123;
	t112 = t122 * t135 + t132 * t179;
	t109 = 0.1e1 / t112 ^ 2;
	t177 = t135 * t123;
	t180 = t133 * t122;
	t111 = t132 * t180 - t177;
	t186 = t111 * t123;
	t108 = 0.1e1 / t112;
	t188 = t108 * t122;
	t146 = t109 * t186 - t188;
	t107 = t111 ^ 2;
	t100 = t107 * t109 + 0.1e1;
	t97 = 0.1e1 / t100;
	t197 = t146 * t97;
	t134 = cos(qJ(3));
	t176 = t135 * t134;
	t117 = atan2(-t176, t132);
	t115 = sin(t117);
	t116 = cos(t117);
	t104 = -t115 * t176 + t116 * t132;
	t101 = 0.1e1 / t104;
	t125 = 0.1e1 / t132;
	t102 = 0.1e1 / t104 ^ 2;
	t126 = 0.1e1 / t132 ^ 2;
	t131 = t135 ^ 2;
	t130 = t134 ^ 2;
	t183 = t126 * t130;
	t120 = t131 * t183 + 0.1e1;
	t118 = 0.1e1 / t120;
	t196 = t118 - 0.1e1;
	t128 = t133 ^ 2;
	t173 = qJD(1) * t135;
	t150 = t130 * t133 * t173;
	t171 = qJD(3) * t134;
	t182 = t128 * t130;
	t170 = qJD(3) * t135;
	t160 = t126 * t170;
	t174 = qJD(1) * t134;
	t162 = t133 * t174;
	t94 = ((t132 * t170 + t162) * t125 + t130 * t160) * t118;
	t155 = -t94 + t170;
	t156 = -t135 * t94 + qJD(3);
	t185 = t116 * t134;
	t88 = t156 * t185 + (t155 * t132 + t162) * t115;
	t192 = t101 * t102 * t88;
	t99 = t102 * t182 + 0.1e1;
	t195 = (-t182 * t192 + (-t128 * t132 * t171 + t150) * t102) / t99 ^ 2;
	t187 = t109 * t111;
	t153 = qJD(1) * t132 + qJD(5);
	t144 = t133 * t171 + t153 * t135;
	t154 = qJD(5) * t132 + qJD(1);
	t148 = t122 * t154;
	t93 = t144 * t123 - t133 * t148;
	t191 = t108 * t109 * t93;
	t147 = t123 * t154;
	t92 = t144 * t122 + t133 * t147;
	t194 = (-t107 * t191 + t92 * t187) / t100 ^ 2;
	t95 = 0.1e1 / t99;
	t193 = t102 * t95;
	t129 = t134 * t130;
	t145 = qJD(3) * (-t126 * t129 - t134) * t125;
	t189 = (-t126 * t150 + t131 * t145) / t120 ^ 2;
	t184 = t125 * t130;
	t181 = t132 * t135;
	t175 = qJD(1) * t133;
	t172 = qJD(3) * t132;
	t169 = -0.2e1 * t195;
	t168 = 0.2e1 * t194;
	t167 = -0.2e1 * t192;
	t166 = t101 * t195;
	t165 = t95 * t172;
	t164 = t134 * t189;
	t163 = t118 * t184;
	t161 = t134 * t173;
	t159 = 0.1e1 + t183;
	t158 = 0.2e1 * t111 * t191;
	t157 = t189 * t198;
	t152 = t135 * t163;
	t151 = t196 * t134 * t115;
	t149 = t159 * t133;
	t143 = -t153 * t133 + t134 * t170;
	t114 = t132 * t177 - t180;
	t113 = t122 * t181 + t179;
	t106 = t159 * t135 * t118;
	t91 = (-t116 * t152 - t151) * t133;
	t90 = t115 * t181 + t185 + (-t115 * t132 - t116 * t176) * t106;
	t89 = -t159 * t157 + (-qJD(1) * t149 + t145 * t198) * t118;
	t1 = [-0.2e1 * t133 * t125 * t164 + (-qJD(3) * t149 + t125 * t161) * t118, 0, t89, 0, 0, 0; (t101 * t165 + (0.2e1 * t166 + (qJD(1) * t91 + t88) * t193) * t134) * t135 + (t91 * t134 * t95 * t167 + (-t91 * t165 + (t91 * t169 + ((t94 * t152 + t196 * t172 + 0.2e1 * t164) * t115 + (t157 * t184 + t94 * t134 + (t129 * t160 + (-t94 + 0.2e1 * t170) * t134) * t118) * t116) * t95 * t133) * t134) * t102 + (t101 + ((t128 - t131) * t116 * t163 - t135 * t151) * t102) * t95 * t174) * t133, 0, (t101 * t95 * t173 + (-0.2e1 * t166 + (-qJD(3) * t90 - t88) * t193) * t133) * t132 + (((qJD(3) * t101 + t90 * t167) * t133 + (t90 * t173 + ((t106 * t175 - t135 * t89) * t116 + ((t106 * t135 - 0.1e1) * t94 + (-t106 + t135) * qJD(3)) * t115) * t133 * t134) * t102) * t95 + (t90 * t169 + ((-t89 - t175) * t115 + (t155 * t106 - t156) * t116) * t95 * t132) * t102 * t133) * t134, 0, 0, 0; (-t108 * t113 + t114 * t187) * t168 + (t114 * t158 + t135 * t108 * t147 + t143 * t188 + (t135 * t111 * t148 - t113 * t93 - t114 * t92 - t143 * t186) * t109) * t97, 0, -t161 * t197 + (t172 * t197 + (t146 * t168 + ((qJD(5) * t108 + t158) * t123 + (-t123 * t92 + (qJD(5) * t111 - t93) * t122) * t109) * t97) * t134) * t133, 0, -0.2e1 * t194 + 0.2e1 * (t109 * t92 * t97 + (-t109 * t194 - t97 * t191) * t111) * t111, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:42:48
	% EndTime: 2019-10-10 00:42:49
	% DurationCPUTime: 1.46s
	% Computational Cost: add. (4821->119), mult. (6168->264), div. (1114->15), fcn. (7752->9), ass. (0->112)
	t157 = sin(qJ(1));
	t158 = cos(qJ(3));
	t207 = t157 * t158;
	t150 = pkin(9) + qJ(5);
	t148 = sin(t150);
	t156 = sin(qJ(3));
	t228 = cos(qJ(1));
	t189 = t228 * t156;
	t149 = cos(t150);
	t208 = t157 * t149;
	t137 = t148 * t189 + t208;
	t206 = t158 * t148;
	t125 = atan2(t137, t206);
	t122 = cos(t125);
	t121 = sin(t125);
	t219 = t121 * t137;
	t116 = t122 * t206 + t219;
	t113 = 0.1e1 / t116;
	t136 = t228 * t148 + t156 * t208;
	t131 = 0.1e1 / t136;
	t145 = 0.1e1 / t148;
	t153 = 0.1e1 / t158;
	t114 = 0.1e1 / t116 ^ 2;
	t132 = 0.1e1 / t136 ^ 2;
	t146 = 0.1e1 / t148 ^ 2;
	t209 = t157 * t148;
	t135 = -t228 * t149 + t156 * t209;
	t130 = t135 ^ 2;
	t111 = t114 * t130 + 0.1e1;
	t184 = qJD(1) * t228;
	t202 = qJD(3) * t158;
	t170 = -t156 * t184 - t157 * t202;
	t167 = t228 * qJD(5) - t170;
	t179 = qJD(5) * t156 + qJD(1);
	t119 = t167 * t148 + t179 * t208;
	t222 = t119 * t114;
	t134 = t137 ^ 2;
	t154 = 0.1e1 / t158 ^ 2;
	t212 = t146 * t154;
	t126 = t134 * t212 + 0.1e1;
	t123 = 0.1e1 / t126;
	t200 = qJD(5) * t158;
	t204 = qJD(3) * t156;
	t171 = -t148 * t204 + t149 * t200;
	t191 = t137 * t212;
	t188 = t228 * t158;
	t176 = qJD(3) * t188;
	t177 = t149 * t189;
	t117 = -qJD(5) * t177 - t148 * t176 - t149 * t184 + (qJD(1) * t156 + qJD(5)) * t209;
	t214 = t145 * t153;
	t194 = t117 * t214;
	t105 = (-t171 * t191 - t194) * t123;
	t169 = -t105 * t137 - t171;
	t101 = (-t105 * t206 - t117) * t121 - t169 * t122;
	t115 = t113 * t114;
	t226 = t101 * t115;
	t227 = 0.1e1 / t111 ^ 2 * (-t130 * t226 + t135 * t222);
	t147 = t145 * t146;
	t152 = t158 ^ 2;
	t155 = t153 / t152;
	t201 = qJD(5) * t149;
	t186 = t154 * t201;
	t225 = (-t117 * t191 + (t146 * t155 * t204 - t147 * t186) * t134) / t126 ^ 2;
	t151 = t157 ^ 2;
	t211 = t151 * t152;
	t193 = t132 * t211;
	t129 = 0.1e1 + t193;
	t168 = -t151 * t156 * t202 + t152 * t157 * t184;
	t120 = t167 * t149 - t179 * t209;
	t221 = t120 * t131 * t132;
	t178 = t211 * t221;
	t224 = (t168 * t132 - t178) / t129 ^ 2;
	t223 = t114 * t135;
	t220 = t121 * t135;
	t218 = t121 * t158;
	t217 = t122 * t135;
	t216 = t122 * t137;
	t215 = t122 * t156;
	t213 = t146 * t149;
	t210 = t156 * t157;
	t205 = qJD(1) * t157;
	t203 = qJD(3) * t157;
	t199 = 0.2e1 * t227;
	t198 = -0.2e1 * t225;
	t197 = 0.2e1 * t224;
	t196 = 0.2e1 * t115 * t135;
	t195 = t114 * t220;
	t192 = t137 * t214;
	t190 = t145 * t154 * t156;
	t187 = t154 * t204;
	t173 = t137 * t190 + t228;
	t112 = t173 * t123;
	t185 = t228 - t112;
	t183 = -0.2e1 * t113 * t227;
	t182 = t114 * t199;
	t181 = 0.2e1 * t153 * t225;
	t180 = -0.2e1 * t135 * t207;
	t175 = t145 * t181;
	t138 = t177 - t209;
	t174 = t137 * t213 - t138 * t145;
	t172 = t132 * t138 * t157 - t228 * t131;
	t166 = t119 * t214 - (t146 * t153 * t201 - t145 * t187) * t135;
	t127 = 0.1e1 / t129;
	t118 = t136 * qJD(1) + qJD(5) * t137 - t149 * t176;
	t109 = 0.1e1 / t111;
	t108 = t174 * t153 * t123;
	t104 = (-t121 + (-t122 * t192 + t121) * t123) * t135;
	t103 = t112 * t216 + (t185 * t218 - t215) * t148;
	t102 = t122 * t149 * t158 + t121 * t138 - (-t121 * t206 + t216) * t108;
	t100 = t173 * t198 + (-t117 * t190 - t205 + (-t146 * t156 * t186 + (0.2e1 * t155 * t156 ^ 2 + t153) * t145 * qJD(3)) * t137) * t123;
	t98 = t174 * t181 + (-t174 * t187 + (t117 * t213 - t118 * t145 + (-t138 * t213 + (0.2e1 * t147 * t149 ^ 2 + t145) * t137) * qJD(5)) * t153) * t123;
	t1 = [-t166 * t123 + t135 * t175, 0, t100, 0, t98, 0; t137 * t183 + (-t117 * t113 + (-t101 * t137 - t104 * t119) * t114) * t109 + (t104 * t182 + (0.2e1 * t104 * t226 + (-t119 * t123 + t119 - (t105 * t123 * t192 + t198) * t135) * t114 * t121 + (-(t137 * t175 - t105) * t223 + (-(t105 + t194) * t135 + t166 * t137) * t114 * t123) * t122) * t109) * t135, 0, t103 * t135 * t182 + (-(t100 * t216 + (-t105 * t219 - t117 * t122) * t112) * t223 + (t101 * t196 - t222) * t103 + (t113 * t207 - (-t112 * t218 + t121 * t188 - t215) * t223) * t201) * t109 + (t183 * t207 + ((-t113 * t203 - (-t185 * qJD(3) + t105) * t195) * t156 + (t113 * t184 + (-t157 * t101 - (-t100 - t205) * t220 - (t185 * t105 - qJD(3)) * t217) * t114) * t158) * t109) * t148, 0, (t102 * t223 - t113 * t136) * t199 + (-t102 * t222 + t120 * t113 + (t102 * t196 - t114 * t136) * t101 - (-t149 * t204 - t148 * t200 + t108 * t117 + t137 * t98 + (t108 * t206 + t138) * t105) * t114 * t217 - (-t118 + (-t105 * t149 - t148 * t98) * t158 - t169 * t108) * t195) * t109, 0; t172 * t158 * t197 + (t172 * t204 + ((-qJD(1) * t131 + 0.2e1 * t138 * t221) * t157 + (t118 * t157 - t228 * t120 - t138 * t184) * t132) * t158) * t127, 0, (t131 * t210 + t149 * t193) * t197 + (0.2e1 * t149 * t178 + t170 * t131 + (qJD(5) * t148 * t211 + t120 * t210 - 0.2e1 * t149 * t168) * t132) * t127, 0, t132 * t180 * t224 + (t180 * t221 + (t119 * t207 + (-t156 * t203 + t158 * t184) * t135) * t132) * t127, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end