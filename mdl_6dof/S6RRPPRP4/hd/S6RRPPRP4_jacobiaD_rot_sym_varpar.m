% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRP4
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
%   Wie in S6RRPPRP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:32
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRP4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRP4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:32:13
	% EndTime: 2019-10-10 09:32:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:32:13
	% EndTime: 2019-10-10 09:32:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:32:13
	% EndTime: 2019-10-10 09:32:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:32:14
	% EndTime: 2019-10-10 09:32:15
	% DurationCPUTime: 1.02s
	% Computational Cost: add. (892->82), mult. (2191->191), div. (456->12), fcn. (2616->9), ass. (0->85)
	t100 = sin(qJ(2));
	t92 = t100 ^ 2;
	t102 = cos(qJ(2));
	t95 = 0.1e1 / t102 ^ 2;
	t141 = t92 * t95;
	t101 = sin(qJ(1));
	t122 = 0.1e1 + t141;
	t93 = t101 ^ 2;
	t90 = t93 * t141 + 0.1e1;
	t88 = 0.1e1 / t90;
	t113 = t122 * t88;
	t74 = t101 * t113;
	t157 = t101 * t74 - 0.1e1;
	t103 = cos(qJ(1));
	t127 = qJD(2) * t103;
	t118 = t100 * t127;
	t134 = t101 * t102;
	t98 = sin(pkin(9));
	t99 = cos(pkin(9));
	t82 = t103 * t98 - t99 * t134;
	t76 = t82 * qJD(1) - t99 * t118;
	t133 = t102 * t103;
	t84 = t101 * t98 + t99 * t133;
	t78 = 0.1e1 / t84;
	t79 = 0.1e1 / t84 ^ 2;
	t80 = t78 * t79;
	t145 = t76 * t80;
	t81 = -t103 * t99 - t98 * t134;
	t75 = t81 * qJD(1) - t98 * t118;
	t146 = t75 * t79;
	t83 = -t101 * t99 + t98 * t133;
	t77 = t83 ^ 2;
	t72 = t77 * t79 + 0.1e1;
	t156 = (-t77 * t145 + t83 * t146) / t72 ^ 2;
	t155 = t100 * t141;
	t143 = t83 * t99;
	t112 = t79 * t143 - t78 * t98;
	t70 = 0.1e1 / t72;
	t154 = t112 * t70;
	t135 = t101 * t100;
	t87 = atan2(-t135, -t102);
	t85 = sin(t87);
	t123 = t85 * t135;
	t86 = cos(t87);
	t69 = -t102 * t86 - t123;
	t66 = 0.1e1 / t69;
	t94 = 0.1e1 / t102;
	t67 = 0.1e1 / t69 ^ 2;
	t153 = 0.2e1 * t100;
	t152 = t88 - 0.1e1;
	t130 = qJD(1) * t103;
	t114 = t101 * t92 * t130;
	t128 = qJD(2) * t102;
	t97 = t103 ^ 2;
	t140 = t92 * t97;
	t129 = qJD(2) * t101;
	t138 = t102 * t85;
	t62 = (-(-t100 * t130 - t101 * t128) * t94 + t129 * t141) * t88;
	t57 = (t62 - t129) * t138 + (-t85 * t130 + (-t101 * t62 + qJD(2)) * t86) * t100;
	t150 = t57 * t66 * t67;
	t65 = t67 * t140 + 0.1e1;
	t151 = (-t140 * t150 + (t100 * t97 * t128 - t114) * t67) / t65 ^ 2;
	t63 = 0.1e1 / t65;
	t148 = t63 * t67;
	t111 = qJD(2) * (t100 + t155) * t94;
	t147 = (t93 * t111 + t95 * t114) / t90 ^ 2;
	t144 = t82 * t83;
	t142 = t88 * t94;
	t137 = t103 * t67;
	t136 = qJD(2) * t74;
	t132 = qJD(1) * t100;
	t131 = qJD(1) * t101;
	t126 = 0.2e1 * t150;
	t125 = t66 * t151;
	t124 = t101 * t142;
	t121 = t100 * t152;
	t120 = t63 * t128;
	t119 = t100 * t129;
	t117 = 0.2e1 * t67 * t151;
	t116 = -0.2e1 * t94 * t147;
	t115 = t92 * t124;
	t61 = (-t86 * t115 + t85 * t121) * t103;
	t59 = (-t101 + t74) * t138 - t157 * t86 * t100;
	t58 = t113 * t130 + 0.2e1 * (t111 * t88 - t122 * t147) * t101;
	t1 = [-t124 * t132 + (qJD(2) * t113 + t100 * t116) * t103, t58, 0, 0, 0, 0; (-t66 * t120 + (0.2e1 * t125 + (qJD(1) * t61 + t57) * t148) * t100) * t101 + (-t61 * t67 * t120 + (t61 * t117 + (t61 * t126 + ((-t62 * t115 - t152 * t128 + t147 * t153) * t85 + (-t62 * t121 + (t92 * t116 + (t153 + t155) * t88 * qJD(2)) * t101) * t86) * t137) * t63) * t100 + (-t66 + (-(t93 - t97) * t92 * t86 * t142 + t152 * t123) * t67) * t63 * t132) * t103, (-t66 * t63 * t131 + (-0.2e1 * t125 + (-qJD(2) * t59 - t57) * t148) * t103) * t102 + (t59 * t103 * t117 + (-t66 * t127 - ((-t101 * t58 - t130 * t74) * t86 + (t157 * t62 + t129 - t136) * t85) * t100 * t137 + (t103 * t126 + t67 * t131) * t59) * t63 - ((t58 - t130) * t85 + (t62 * t74 + qJD(2) + (-t62 - t136) * t101) * t86) * t133 * t148) * t100, 0, 0, 0, 0; 0.2e1 * (t79 * t144 - t78 * t81) * t156 + ((-t83 * qJD(1) + t98 * t119) * t78 + 0.2e1 * t144 * t145 + (-t81 * t76 - (-t84 * qJD(1) + t99 * t119) * t83 - t82 * t75) * t79) * t70, t102 * t127 * t154 + (-t131 * t154 + (-0.2e1 * t112 * t156 + (t99 * t146 + (-0.2e1 * t80 * t143 + t79 * t98) * t76) * t70) * t103) * t100, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:32:14
	% EndTime: 2019-10-10 09:32:14
	% DurationCPUTime: 0.80s
	% Computational Cost: add. (926->76), mult. (3228->188), div. (613->15), fcn. (4191->9), ass. (0->85)
	t120 = cos(pkin(9));
	t172 = 0.2e1 * t120;
	t121 = sin(qJ(2));
	t115 = 0.1e1 / t121;
	t119 = sin(pkin(9));
	t124 = cos(qJ(1));
	t149 = t124 * t120;
	t122 = sin(qJ(1));
	t123 = cos(qJ(2));
	t152 = t122 * t123;
	t104 = t119 * t152 + t149;
	t112 = 0.1e1 / t119;
	t158 = t104 * t112;
	t140 = t115 * t158;
	t154 = t121 * t119;
	t94 = atan2(-t104, t154);
	t90 = sin(t94);
	t91 = cos(t94);
	t113 = 0.1e1 / t119 ^ 2;
	t116 = 0.1e1 / t121 ^ 2;
	t99 = t104 ^ 2;
	t97 = t99 * t116 * t113 + 0.1e1;
	t92 = 0.1e1 / t97;
	t171 = (t91 * t140 + t90) * t92 - t90;
	t86 = -t90 * t104 + t91 * t154;
	t83 = 0.1e1 / t86;
	t150 = t124 * t119;
	t138 = t123 * t150;
	t107 = -t122 * t120 + t138;
	t75 = t171 * t107;
	t170 = 0.2e1 * t75;
	t153 = t122 * t119;
	t108 = t123 * t149 + t153;
	t101 = 0.1e1 / t108;
	t102 = 0.1e1 / t108 ^ 2;
	t84 = 0.1e1 / t86 ^ 2;
	t100 = t107 ^ 2;
	t145 = qJD(2) * t124;
	t136 = t121 * t145;
	t87 = t104 * qJD(1) + t119 * t136;
	t166 = t84 * t87;
	t146 = qJD(2) * t123;
	t161 = t121 * t90;
	t164 = t104 * t91;
	t137 = t116 * t146;
	t148 = qJD(1) * t122;
	t89 = -qJD(2) * t121 * t153 + qJD(1) * t138 - t120 * t148;
	t76 = (t104 * t137 - t115 * t89) * t92 * t112;
	t73 = -t76 * t164 - t90 * t89 + (t91 * t146 - t76 * t161) * t119;
	t168 = t73 * t83 * t84;
	t80 = t100 * t84 + 0.1e1;
	t169 = (-t100 * t168 - t107 * t166) / t80 ^ 2;
	t114 = t121 ^ 2;
	t117 = t115 / t114;
	t167 = (t104 * t116 * t89 - t117 * t99 * t146) * t113 / t97 ^ 2;
	t106 = -t120 * t152 + t150;
	t133 = t120 * t136;
	t88 = t106 * qJD(1) - t133;
	t165 = t101 * t102 * t88;
	t163 = t107 * t91;
	t162 = t115 * t92;
	t160 = t90 * t107;
	t155 = t116 * t123;
	t82 = (t155 * t158 + t122) * t92;
	t159 = t122 - t82;
	t157 = t106 * t124;
	t118 = t124 ^ 2;
	t156 = t114 * t118;
	t151 = t123 * t124;
	t147 = qJD(1) * t124;
	t130 = t114 * t122 * t147 - t118 * t121 * t146;
	t134 = t156 * t165;
	t139 = t102 * t156;
	t98 = 0.1e1 + t139;
	t144 = 0.2e1 * (-t130 * t102 - t134) / t98 ^ 2;
	t143 = -0.2e1 * t167;
	t142 = t84 * t169;
	t141 = t84 * t160;
	t135 = 0.2e1 * t83 * t169;
	t131 = 0.2e1 * t115 * t167 + t92 * t137;
	t95 = 0.1e1 / t98;
	t78 = 0.1e1 / t80;
	t74 = -t82 * t164 + (t123 * t91 + t159 * t161) * t119;
	t72 = t122 * t143 + t92 * t147 + (t89 * t92 * t155 + (t143 * t155 + (-0.2e1 * t117 * t123 ^ 2 - t115) * t92 * qJD(2)) * t104) * t112;
	t1 = [(t131 * t107 + t87 * t162) * t112, t72, 0, 0, 0, 0; t104 * t135 + (-t89 * t83 + (t104 * t73 + t75 * t87) * t84) * t78 + (t142 * t170 + (t168 * t170 - (-t76 * t92 * t140 + t143) * t141 - ((t92 - 0.1e1) * t76 + (-t131 * t104 + t89 * t162) * t112) * t84 * t163 + t171 * t166) * t78) * t107, t74 * t78 * t166 + (-(-t82 * t91 * t89 + (t76 * t82 * t90 - t72 * t91) * t104) * t84 * t78 + 0.2e1 * (t78 * t168 + t142) * t74) * t107 + (t124 * t135 * t121 + ((-t83 * t145 - (t159 * qJD(2) - t76) * t141) * t123 + (t83 * t148 + (t124 * t73 - (-t72 + t147) * t160 - (t159 * t76 - qJD(2)) * t163) * t84) * t121) * t78) * t119, 0, 0, 0, 0; (-t101 * t122 - t102 * t157) * t121 * t144 + (-0.2e1 * t121 * t157 * t165 + (t121 * t147 + t122 * t146) * t101 + (((-t88 + t133) * t122 - t108 * t147) * t121 + (-t121 * t148 + t123 * t145) * t106) * t102) * t95, (t101 * t151 + t120 * t139) * t144 + (t134 * t172 + (t123 * t148 + t136) * t101 + (t130 * t172 + t88 * t151) * t102) * t95, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:32:14
	% EndTime: 2019-10-10 09:32:15
	% DurationCPUTime: 1.28s
	% Computational Cost: add. (926->105), mult. (3345->234), div. (468->12), fcn. (4023->11), ass. (0->103)
	t155 = sin(qJ(2));
	t146 = t155 ^ 2;
	t158 = cos(qJ(2));
	t149 = 0.1e1 / t158 ^ 2;
	t204 = t146 * t149;
	t156 = sin(qJ(1));
	t147 = t156 ^ 2;
	t144 = t147 * t204 + 0.1e1;
	t148 = 0.1e1 / t158;
	t201 = t148 * t155;
	t222 = t155 * t204;
	t167 = qJD(2) * (t148 * t222 + t201);
	t159 = cos(qJ(1));
	t195 = qJD(1) * t159;
	t202 = t146 * t156;
	t173 = t195 * t202;
	t207 = (t147 * t167 + t149 * t173) / t144 ^ 2;
	t223 = -0.2e1 * t207;
	t177 = 0.1e1 + t204;
	t221 = t156 * t177;
	t154 = sin(qJ(5));
	t157 = cos(qJ(5));
	t191 = qJD(5) * t159;
	t196 = qJD(1) * t156;
	t220 = t154 * t196 - t157 * t191;
	t219 = t154 * t191 + t157 * t196;
	t152 = sin(pkin(9));
	t153 = cos(pkin(9));
	t197 = t158 * t159;
	t137 = t152 * t197 - t156 * t153;
	t138 = t156 * t152 + t153 * t197;
	t121 = t137 * t154 + t138 * t157;
	t115 = 0.1e1 / t121;
	t199 = t156 * t155;
	t143 = atan2(t199, t158);
	t140 = cos(t143);
	t139 = sin(t143);
	t186 = t139 * t199;
	t125 = t140 * t158 + t186;
	t122 = 0.1e1 / t125;
	t116 = 0.1e1 / t121 ^ 2;
	t123 = 0.1e1 / t125 ^ 2;
	t218 = 0.2e1 * t155;
	t141 = 0.1e1 / t144;
	t217 = t141 - 0.1e1;
	t198 = t156 * t158;
	t135 = -t152 * t198 - t153 * t159;
	t192 = qJD(2) * t159;
	t180 = t155 * t192;
	t128 = t135 * qJD(1) - t152 * t180;
	t136 = t152 * t159 - t153 * t198;
	t129 = t136 * qJD(1) - t153 * t180;
	t104 = t121 * qJD(5) - t128 * t157 + t129 * t154;
	t170 = t137 * t157 - t138 * t154;
	t114 = t170 ^ 2;
	t108 = t114 * t116 + 0.1e1;
	t210 = t116 * t170;
	t105 = t170 * qJD(5) + t128 * t154 + t129 * t157;
	t117 = t115 * t116;
	t213 = t105 * t117;
	t216 = 0.1e1 / t108 ^ 2 * (-t104 * t210 - t114 * t213);
	t151 = t159 ^ 2;
	t203 = t146 * t151;
	t113 = t123 * t203 + 0.1e1;
	t193 = qJD(2) * t158;
	t183 = t155 * t195;
	t194 = qJD(2) * t156;
	t110 = ((t156 * t193 + t183) * t148 + t194 * t204) * t141;
	t205 = t140 * t155;
	t101 = (t110 * t156 - qJD(2)) * t205 + (t183 + (-t110 + t194) * t158) * t139;
	t214 = t101 * t122 * t123;
	t215 = (-t203 * t214 + (t151 * t155 * t193 - t173) * t123) / t113 ^ 2;
	t212 = t110 * t139;
	t211 = t110 * t155;
	t168 = -t152 * t154 - t153 * t157;
	t200 = t155 * t159;
	t133 = t168 * t200;
	t209 = t116 * t133;
	t208 = t123 * t155;
	t127 = t141 * t221;
	t206 = t127 * t156;
	t190 = 0.2e1 * t216;
	t189 = -0.2e1 * t214;
	t188 = -0.2e1 * t117 * t170;
	t187 = t123 * t200;
	t185 = t141 * t146 * t148;
	t181 = t155 * t194;
	t176 = -0.2e1 * t155 * t215;
	t175 = t148 * t223;
	t174 = t156 * t185;
	t172 = t177 * t159;
	t171 = t135 * t157 - t136 * t154;
	t119 = t135 * t154 + t136 * t157;
	t169 = t152 * t157 - t153 * t154;
	t132 = t169 * t200;
	t131 = -t138 * qJD(1) + t153 * t181;
	t130 = -t137 * qJD(1) + t152 * t181;
	t111 = 0.1e1 / t113;
	t109 = (-t217 * t155 * t139 + t140 * t174) * t159;
	t106 = 0.1e1 / t108;
	t103 = t139 * t198 - t205 + (-t139 * t158 + t140 * t199) * t127;
	t102 = t221 * t223 + (qJD(1) * t172 + 0.2e1 * t156 * t167) * t141;
	t1 = [t175 * t200 + (qJD(2) * t172 - t196 * t201) * t141, t102, 0, 0, 0, 0; (t122 * t176 + (t122 * t193 + (-qJD(1) * t109 - t101) * t208) * t111) * t156 + (t123 * t176 * t109 + (((-t110 * t174 - t217 * t193 + t207 * t218) * t139 + (t175 * t202 + t211 + (-t211 + (t218 + t222) * t194) * t141) * t140) * t187 + (t123 * t193 + t155 * t189) * t109 + (t122 + ((-t147 + t151) * t140 * t185 + t217 * t186) * t123) * t155 * qJD(1)) * t111) * t159, 0.2e1 * (-t103 * t208 + t122 * t158) * t159 * t215 + ((t122 * t196 + (qJD(2) * t103 + t101) * t159 * t123) * t158 + (t122 * t192 + (t102 * t140 * t156 - t139 * t194 - t206 * t212 + t212 + (qJD(2) * t139 + t140 * t195) * t127) * t187 + (-t123 * t196 + t159 * t189) * t103 + ((-t102 + t195) * t139 + ((-0.1e1 + t206) * qJD(2) + (-t127 + t156) * t110) * t140) * t123 * t197) * t155) * t111, 0, 0, 0, 0; (t115 * t171 - t119 * t210) * t190 + ((t119 * qJD(5) - t130 * t157 + t131 * t154) * t115 + t119 * t105 * t188 + (t171 * t105 + (t171 * qJD(5) + t130 * t154 + t131 * t157) * t170 - t119 * t104) * t116) * t106, (-t115 * t132 - t170 * t209) * t190 + (-t104 * t209 + (-t116 * t132 + t133 * t188) * t105 + (t169 * t115 + t168 * t210) * t158 * t192 + ((t220 * t115 + t219 * t210) * t153 + (-t219 * t115 + t220 * t210) * t152) * t155) * t106, 0, 0, -0.2e1 * t216 - 0.2e1 * (t104 * t106 * t116 - (-t106 * t213 - t116 * t216) * t170) * t170, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:32:14
	% EndTime: 2019-10-10 09:32:16
	% DurationCPUTime: 2.18s
	% Computational Cost: add. (4632->128), mult. (14196->269), div. (690->12), fcn. (17956->11), ass. (0->121)
	t217 = cos(pkin(9));
	t216 = sin(pkin(9));
	t223 = cos(qJ(1));
	t267 = t223 * t216;
	t220 = sin(qJ(1));
	t222 = cos(qJ(2));
	t269 = t220 * t222;
	t207 = t217 * t269 - t267;
	t218 = sin(qJ(5));
	t221 = cos(qJ(5));
	t240 = t216 * t269 + t217 * t223;
	t187 = t207 * t221 + t240 * t218;
	t268 = t222 * t223;
	t271 = t220 * t216;
	t208 = t217 * t268 + t271;
	t219 = sin(qJ(2));
	t273 = t219 * t220;
	t256 = qJD(2) * t273;
	t195 = t208 * qJD(1) - t217 * t256;
	t257 = t222 * t267;
	t266 = qJD(1) * t220;
	t237 = qJD(1) * t257 - t216 * t256 - t217 * t266;
	t165 = t187 * qJD(5) + t195 * t218 - t237 * t221;
	t296 = t207 * t218 - t240 * t221;
	t303 = qJD(5) * t296 - t195 * t221 - t237 * t218;
	t270 = t220 * t217;
	t241 = t257 - t270;
	t190 = t208 * t218 - t241 * t221;
	t302 = -0.2e1 * t190;
	t191 = t208 * t221 + t241 * t218;
	t183 = 0.1e1 / t191 ^ 2;
	t214 = t219 ^ 2;
	t215 = t223 ^ 2;
	t274 = t214 * t215;
	t176 = t183 * t274 + 0.1e1;
	t264 = qJD(2) * t222;
	t265 = qJD(1) * t223;
	t263 = qJD(2) * t223;
	t255 = t219 * t263;
	t194 = -t207 * qJD(1) - t217 * t255;
	t230 = -t240 * qJD(1) - t216 * t255;
	t164 = -t190 * qJD(5) + t194 * t221 + t230 * t218;
	t182 = 0.1e1 / t191;
	t288 = t164 * t182 * t183;
	t301 = (-t274 * t288 + (-t214 * t220 * t265 + t215 * t219 * t264) * t183) / t176 ^ 2;
	t247 = t216 * t221 - t217 * t218;
	t206 = t247 * t222;
	t232 = qJD(2) * t206;
	t246 = t216 * t218 + t217 * t221;
	t234 = qJD(5) * t246;
	t178 = t219 * t234 - t232;
	t299 = t247 * t219;
	t202 = 0.1e1 / t299 ^ 2;
	t300 = t178 * t202;
	t173 = atan2(-t296, -t299);
	t168 = sin(t173);
	t169 = cos(t173);
	t162 = -t168 * t296 - t169 * t299;
	t159 = 0.1e1 / t162;
	t201 = 0.1e1 / t299;
	t160 = 0.1e1 / t162 ^ 2;
	t295 = -0.2e1 * t296;
	t294 = 0.2e1 * t190;
	t181 = t190 ^ 2;
	t158 = t160 * t181 + 0.1e1;
	t163 = t191 * qJD(5) + t194 * t218 - t230 * t221;
	t289 = t160 * t190;
	t180 = t296 ^ 2;
	t172 = t180 * t202 + 0.1e1;
	t170 = 0.1e1 / t172;
	t279 = t296 * t202;
	t245 = t165 * t201 + t178 * t279;
	t151 = t245 * t170;
	t248 = t168 * t299 - t169 * t296;
	t147 = t248 * t151 - t165 * t168 + t169 * t178;
	t292 = t147 * t159 * t160;
	t293 = (t163 * t289 - t181 * t292) / t158 ^ 2;
	t283 = t201 * t300;
	t291 = (t165 * t279 + t180 * t283) / t172 ^ 2;
	t156 = 0.1e1 / t158;
	t290 = t156 * t160;
	t287 = t168 * t190;
	t286 = t169 * t190;
	t174 = 0.1e1 / t176;
	t284 = t174 * t183;
	t282 = t183 * t219;
	t281 = t183 * t222;
	t280 = t296 * t201;
	t272 = t219 * t223;
	t262 = 0.2e1 * t293;
	t261 = -0.2e1 * t291;
	t260 = 0.2e1 * t301;
	t259 = t201 * t291;
	t258 = t183 * t272;
	t254 = t222 * t263;
	t253 = t292 * t294;
	t252 = t283 * t295;
	t251 = t273 * t284;
	t250 = t272 * t288;
	t249 = t258 * t301;
	t238 = t246 * t223;
	t198 = t219 * t238;
	t244 = -t182 * t222 - t198 * t282;
	t205 = t246 * t219;
	t243 = t187 * t201 + t205 * t279;
	t196 = (t218 * t270 - t221 * t271) * t219;
	t242 = -t196 * t201 - t206 * t279;
	t236 = qJD(1) * t247;
	t235 = qJD(2) * t247;
	t233 = -t168 + (-t169 * t280 + t168) * t170;
	t197 = t223 * t299;
	t179 = qJD(5) * t299 + t246 * t264;
	t177 = t219 * t235 + t222 * t234;
	t167 = -t220 * t232 + (t220 * t234 - t223 * t236) * t219;
	t155 = t242 * t170;
	t153 = t243 * t170;
	t149 = t248 * t155 + t168 * t196 - t169 * t206;
	t148 = t248 * t153 - t168 * t187 + t169 * t205;
	t146 = t242 * t261 + (t206 * t252 - t167 * t201 + (-t165 * t206 + t177 * t296 - t178 * t196) * t202) * t170;
	t145 = t243 * t261 + (-t205 * t252 - t303 * t201 + (t165 * t205 + t178 * t187 + t179 * t296) * t202) * t170;
	t1 = [-t259 * t294 + (t163 * t201 + t190 * t300) * t170, t146, 0, 0, t145, 0; -(-t147 * t290 - 0.2e1 * t159 * t293) * t296 + (-t165 * t159 - (t233 * t163 + ((t151 * t170 * t280 + t261) * t168 + (-t259 * t295 - t151 + (t151 - t245) * t170) * t169) * t190) * t289) * t156 + (t156 * t253 - t163 * t290 + t289 * t262) * t233 * t190, (t149 * t289 - t159 * t197) * t262 + (t149 * t253 + (-t197 * t147 - t149 * t163 - (-t146 * t296 - t155 * t165 + t177 + (t155 * t299 + t196) * t151) * t286 - (t146 * t299 - t155 * t178 + t167 + (t155 * t296 + t206) * t151) * t287) * t160 + (t235 * t268 + (-qJD(5) * t238 - t220 * t236) * t219) * t159) * t156, 0, 0, (t148 * t289 - t159 * t191) * t262 + (t148 * t253 + t164 * t159 + (-t191 * t147 - t148 * t163 - (-t145 * t296 - t153 * t165 + t179 + (t153 * t299 - t187) * t151) * t286 - (t145 * t299 - t153 * t178 + t303 + (t153 * t296 - t205) * t151) * t287) * t160) * t156, 0; t164 * t251 + t182 * t260 * t273 - (qJD(1) * t251 - t254 * t284 + 0.2e1 * t249) * t187 + (-t303 * t258 + (-t219 * t265 - t220 * t264) * t182 - 0.2e1 * t187 * t250) * t174, t244 * t223 * t260 + (t244 * t266 + ((qJD(2) * t198 - t164) * t281 + (-0.2e1 * t198 * t288 + (t247 * t223 * qJD(5) - t246 * t266) * t282 + (t238 * t281 - t182) * qJD(2)) * t219) * t223) * t174, 0, 0, t249 * t302 + (t250 * t302 + (t163 * t272 + (-t219 * t266 + t254) * t190) * t183) * t174, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end