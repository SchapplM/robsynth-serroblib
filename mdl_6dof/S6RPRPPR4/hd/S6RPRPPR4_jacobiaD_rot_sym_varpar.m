% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR4
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
%   Wie in S6RPRPPR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:20
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:29
	% DurationCPUTime: 1.00s
	% Computational Cost: add. (2086->84), mult. (2191->187), div. (456->12), fcn. (2616->9), ass. (0->88)
	t104 = pkin(9) + qJ(3);
	t103 = cos(t104);
	t100 = 0.1e1 / t103 ^ 2;
	t102 = sin(t104);
	t98 = t102 ^ 2;
	t151 = t100 * t98;
	t109 = sin(qJ(1));
	t127 = 0.1e1 + t151;
	t105 = t109 ^ 2;
	t96 = t105 * t151 + 0.1e1;
	t94 = 0.1e1 / t96;
	t120 = t127 * t94;
	t77 = t109 * t120;
	t168 = t109 * t77 - 0.1e1;
	t108 = cos(pkin(10));
	t110 = cos(qJ(1));
	t136 = qJD(3) * t110;
	t125 = t102 * t136;
	t142 = t109 * t108;
	t107 = sin(pkin(10));
	t146 = t107 * t110;
	t88 = -t103 * t142 + t146;
	t82 = qJD(1) * t88 - t108 * t125;
	t143 = t109 * t107;
	t145 = t108 * t110;
	t90 = t103 * t145 + t143;
	t85 = 0.1e1 / t90 ^ 2;
	t167 = t82 * t85;
	t166 = t102 * t151;
	t89 = t103 * t146 - t142;
	t153 = t85 * t89;
	t83 = t89 ^ 2;
	t80 = t83 * t85 + 0.1e1;
	t78 = 0.1e1 / t80;
	t84 = 0.1e1 / t90;
	t165 = (-t107 * t84 + t108 * t153) * t78;
	t144 = t109 * t102;
	t93 = atan2(-t144, -t103);
	t91 = sin(t93);
	t130 = t91 * t144;
	t92 = cos(t93);
	t75 = -t103 * t92 - t130;
	t72 = 0.1e1 / t75;
	t99 = 0.1e1 / t103;
	t73 = 0.1e1 / t75 ^ 2;
	t164 = 0.2e1 * t102;
	t163 = t94 - 0.1e1;
	t106 = t110 ^ 2;
	t139 = qJD(1) * t110;
	t121 = t109 * t98 * t139;
	t138 = qJD(3) * t103;
	t128 = t73 * t138;
	t137 = qJD(3) * t109;
	t150 = t103 * t91;
	t68 = (-(-t102 * t139 - t103 * t137) * t99 + t137 * t151) * t94;
	t63 = (t68 - t137) * t150 + (-t91 * t139 + (-t109 * t68 + qJD(3)) * t92) * t102;
	t161 = t63 * t72 * t73;
	t156 = t73 * t98;
	t71 = t106 * t156 + 0.1e1;
	t162 = (-t73 * t121 + (t102 * t128 - t161 * t98) * t106) / t71 ^ 2;
	t154 = t84 * t167;
	t87 = -t103 * t143 - t145;
	t81 = qJD(1) * t87 - t107 * t125;
	t160 = (t153 * t81 - t154 * t83) / t80 ^ 2;
	t69 = 0.1e1 / t71;
	t158 = t69 * t73;
	t157 = t72 * t69;
	t118 = qJD(3) * (t102 + t166) * t99;
	t155 = (t100 * t121 + t105 * t118) / t96 ^ 2;
	t152 = t94 * t99;
	t148 = t110 * t73;
	t147 = qJD(3) * t77;
	t141 = qJD(1) * t102;
	t140 = qJD(1) * t109;
	t135 = 0.2e1 * t161;
	t134 = 0.2e1 * t160;
	t133 = t72 * t162;
	t132 = t89 * t154;
	t131 = t109 * t152;
	t129 = t102 * t163;
	t126 = t102 * t137;
	t124 = 0.2e1 * t73 * t162;
	t123 = -0.2e1 * t99 * t155;
	t122 = t98 * t131;
	t67 = (-t122 * t92 + t129 * t91) * t110;
	t65 = (-t109 + t77) * t150 - t168 * t92 * t102;
	t64 = t120 * t139 + 0.2e1 * (t118 * t94 - t127 * t155) * t109;
	t1 = [-t131 * t141 + (qJD(3) * t120 + t102 * t123) * t110, 0, t64, 0, 0, 0; (-t138 * t157 + (0.2e1 * t133 + (qJD(1) * t67 + t63) * t158) * t102) * t109 + (t67 * t124 * t102 + (-t67 * t128 + (t67 * t135 + ((-t122 * t68 - t138 * t163 + t155 * t164) * t91 + (-t68 * t129 + (t98 * t123 + (t164 + t166) * t94 * qJD(3)) * t109) * t92) * t148) * t102 + (-t72 - (t105 - t106) * t92 * t152 * t156 + t163 * t73 * t130) * t141) * t69) * t110, 0, (-t140 * t157 + (-0.2e1 * t133 + (-qJD(3) * t65 - t63) * t158) * t110) * t103 + (t65 * t110 * t124 + (-t72 * t136 + (t110 * t135 + t140 * t73) * t65 + (-((-t109 * t64 - t139 * t77) * t92 + (t168 * t68 + t137 - t147) * t91) * t102 - ((t64 - t139) * t91 + (t68 * t77 + qJD(3) + (-t68 - t147) * t109) * t92) * t103) * t148) * t69) * t102, 0, 0, 0; (t153 * t88 - t84 * t87) * t134 + ((-qJD(1) * t89 + t107 * t126) * t84 + 0.2e1 * t88 * t132 + (-t87 * t82 - (-qJD(1) * t90 + t108 * t126) * t89 - t88 * t81) * t85) * t78, 0, t103 * t136 * t165 + (-t140 * t165 + ((t134 * t84 + t78 * t167) * t107 + (-0.2e1 * t153 * t160 + (t81 * t85 - 0.2e1 * t132) * t78) * t108) * t110) * t102, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:29
	% DurationCPUTime: 0.78s
	% Computational Cost: add. (2051->78), mult. (3228->189), div. (613->15), fcn. (4191->9), ass. (0->87)
	t130 = cos(pkin(10));
	t181 = 0.2e1 * t130;
	t127 = pkin(9) + qJ(3);
	t123 = sin(t127);
	t120 = 0.1e1 / t123;
	t124 = cos(t127);
	t132 = cos(qJ(1));
	t159 = t132 * t130;
	t129 = sin(pkin(10));
	t131 = sin(qJ(1));
	t162 = t131 * t129;
	t111 = t124 * t162 + t159;
	t125 = 0.1e1 / t129;
	t167 = t111 * t125;
	t149 = t120 * t167;
	t163 = t123 * t129;
	t101 = atan2(-t111, t163);
	t97 = sin(t101);
	t98 = cos(t101);
	t106 = t111 ^ 2;
	t121 = 0.1e1 / t123 ^ 2;
	t126 = 0.1e1 / t129 ^ 2;
	t104 = t106 * t121 * t126 + 0.1e1;
	t99 = 0.1e1 / t104;
	t180 = (t98 * t149 + t97) * t99 - t97;
	t93 = -t97 * t111 + t98 * t163;
	t90 = 0.1e1 / t93;
	t160 = t132 * t129;
	t147 = t124 * t160;
	t161 = t131 * t130;
	t114 = t147 - t161;
	t82 = t180 * t114;
	t179 = 0.2e1 * t82;
	t115 = t124 * t159 + t162;
	t108 = 0.1e1 / t115;
	t109 = 0.1e1 / t115 ^ 2;
	t91 = 0.1e1 / t93 ^ 2;
	t107 = t114 ^ 2;
	t154 = qJD(3) * t132;
	t144 = t123 * t154;
	t94 = t111 * qJD(1) + t129 * t144;
	t176 = t91 * t94;
	t156 = qJD(3) * t124;
	t170 = t123 * t97;
	t173 = t111 * t98;
	t146 = t121 * t156;
	t155 = qJD(3) * t131;
	t158 = qJD(1) * t131;
	t96 = qJD(1) * t147 - t130 * t158 - t155 * t163;
	t83 = (t111 * t146 - t120 * t96) * t99 * t125;
	t80 = -t83 * t173 - t97 * t96 + (t98 * t156 - t83 * t170) * t129;
	t177 = t80 * t90 * t91;
	t87 = t107 * t91 + 0.1e1;
	t178 = (-t107 * t177 - t114 * t176) / t87 ^ 2;
	t119 = t123 ^ 2;
	t122 = t120 / t119;
	t175 = 0.1e1 / t104 ^ 2 * (-t106 * t122 * t156 + t111 * t121 * t96) * t126;
	t113 = -t124 * t161 + t160;
	t140 = t130 * t144;
	t95 = t113 * qJD(1) - t140;
	t174 = t108 * t109 * t95;
	t172 = t114 * t98;
	t171 = t120 * t99;
	t169 = t97 * t114;
	t164 = t121 * t124;
	t89 = (t164 * t167 + t131) * t99;
	t168 = t131 - t89;
	t166 = t113 * t132;
	t128 = t132 ^ 2;
	t165 = t119 * t128;
	t157 = qJD(1) * t132;
	t153 = -0.2e1 * t175;
	t148 = t109 * t165;
	t105 = 0.1e1 + t148;
	t141 = t119 * t131 * t157;
	t142 = t165 * t174;
	t145 = qJD(3) * t123 * t128;
	t152 = 0.2e1 / t105 ^ 2 * (-t142 + (t124 * t145 - t141) * t109);
	t151 = t91 * t178;
	t150 = t91 * t169;
	t143 = 0.2e1 * t90 * t178;
	t138 = 0.2e1 * t120 * t175 + t99 * t146;
	t102 = 0.1e1 / t105;
	t85 = 0.1e1 / t87;
	t81 = -t89 * t173 + (t124 * t98 + t168 * t170) * t129;
	t79 = t131 * t153 + t99 * t157 + (t96 * t99 * t164 + (t153 * t164 + (-0.2e1 * t122 * t124 ^ 2 - t120) * t99 * qJD(3)) * t111) * t125;
	t1 = [(t138 * t114 + t94 * t171) * t125, 0, t79, 0, 0, 0; t111 * t143 + (-t96 * t90 + (t111 * t80 + t82 * t94) * t91) * t85 + (t151 * t179 + (t177 * t179 - (-t83 * t99 * t149 + t153) * t150 - ((t99 - 0.1e1) * t83 + (-t138 * t111 + t96 * t171) * t125) * t91 * t172 + t180 * t176) * t85) * t114, 0, t81 * t85 * t176 + (-(-t89 * t98 * t96 + (t83 * t89 * t97 - t79 * t98) * t111) * t91 * t85 + 0.2e1 * (t85 * t177 + t151) * t81) * t114 + (t132 * t143 * t123 + ((-t90 * t154 - (t168 * qJD(3) - t83) * t150) * t124 + (t90 * t158 + (t132 * t80 - (-t79 + t157) * t169 - (t168 * t83 - qJD(3)) * t172) * t91) * t123) * t85) * t129, 0, 0, 0; (-t108 * t131 - t109 * t166) * t123 * t152 + (-0.2e1 * t123 * t166 * t174 + (t123 * t157 + t124 * t155) * t108 + (((-t95 + t140) * t131 - t115 * t157) * t123 + (-t123 * t158 + t124 * t154) * t113) * t109) * t102, 0, (t108 * t124 * t132 + t130 * t148) * t152 + (t142 * t181 + (t124 * t158 + t144) * t108 + (t141 * t181 + (-0.2e1 * t130 * t145 + t132 * t95) * t124) * t109) * t102, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:29
	% DurationCPUTime: 1.20s
	% Computational Cost: add. (2357->105), mult. (3345->232), div. (468->12), fcn. (4023->11), ass. (0->107)
	t165 = pkin(9) + qJ(3);
	t163 = sin(t165);
	t159 = t163 ^ 2;
	t164 = cos(t165);
	t161 = 0.1e1 / t164 ^ 2;
	t220 = t159 * t161;
	t171 = sin(qJ(1));
	t166 = t171 ^ 2;
	t157 = t166 * t220 + 0.1e1;
	t160 = 0.1e1 / t164;
	t217 = t160 * t163;
	t239 = t163 * t220;
	t181 = qJD(3) * (t160 * t239 + t217);
	t173 = cos(qJ(1));
	t209 = qJD(1) * t173;
	t218 = t159 * t171;
	t187 = t209 * t218;
	t225 = (t161 * t187 + t166 * t181) / t157 ^ 2;
	t240 = -0.2e1 * t225;
	t191 = 0.1e1 + t220;
	t238 = t171 * t191;
	t154 = 0.1e1 / t157;
	t196 = t163 * t209;
	t207 = qJD(3) * t171;
	t123 = ((t164 * t207 + t196) * t160 + t207 * t220) * t154;
	t237 = t123 - t207;
	t170 = sin(qJ(6));
	t172 = cos(qJ(6));
	t205 = qJD(6) * t173;
	t210 = qJD(1) * t171;
	t236 = t170 * t210 - t172 * t205;
	t235 = t170 * t205 + t172 * t210;
	t215 = t171 * t163;
	t156 = atan2(t215, t164);
	t153 = cos(t156);
	t152 = sin(t156);
	t200 = t152 * t215;
	t130 = t153 * t164 + t200;
	t127 = 0.1e1 / t130;
	t168 = sin(pkin(10));
	t212 = t173 * t168;
	t169 = cos(pkin(10));
	t213 = t171 * t169;
	t150 = t164 * t212 - t213;
	t211 = t173 * t169;
	t214 = t171 * t168;
	t151 = t164 * t211 + t214;
	t140 = t150 * t170 + t151 * t172;
	t134 = 0.1e1 / t140;
	t128 = 0.1e1 / t130 ^ 2;
	t135 = 0.1e1 / t140 ^ 2;
	t234 = 0.2e1 * t163;
	t233 = t154 - 0.1e1;
	t167 = t173 ^ 2;
	t219 = t159 * t167;
	t126 = t128 * t219 + 0.1e1;
	t208 = qJD(3) * t164;
	t221 = t153 * t163;
	t114 = (t123 * t171 - qJD(3)) * t221 + (-t237 * t164 + t196) * t152;
	t230 = t114 * t127 * t128;
	t232 = (-t219 * t230 + (t163 * t167 * t208 - t187) * t128) / t126 ^ 2;
	t148 = -t164 * t214 - t211;
	t206 = qJD(3) * t173;
	t194 = t163 * t206;
	t141 = t148 * qJD(1) - t168 * t194;
	t149 = -t164 * t213 + t212;
	t142 = t149 * qJD(1) - t169 * t194;
	t117 = t140 * qJD(6) - t141 * t172 + t142 * t170;
	t184 = t150 * t172 - t151 * t170;
	t133 = t184 ^ 2;
	t122 = t133 * t135 + 0.1e1;
	t224 = t135 * t184;
	t118 = t184 * qJD(6) + t141 * t170 + t142 * t172;
	t136 = t134 * t135;
	t229 = t118 * t136;
	t231 = (-t117 * t224 - t133 * t229) / t122 ^ 2;
	t228 = t123 * t163;
	t227 = t128 * t163;
	t226 = t128 * t173;
	t182 = -t168 * t170 - t169 * t172;
	t216 = t163 * t173;
	t146 = t182 * t216;
	t223 = t135 * t146;
	t222 = t152 * t171;
	t204 = 0.2e1 * t231;
	t203 = -0.2e1 * t230;
	t202 = -0.2e1 * t136 * t184;
	t201 = t128 * t216;
	t199 = t154 * t159 * t160;
	t195 = t163 * t207;
	t190 = -0.2e1 * t163 * t232;
	t189 = t160 * t240;
	t188 = t171 * t199;
	t186 = t191 * t173;
	t185 = t148 * t172 - t149 * t170;
	t138 = t148 * t170 + t149 * t172;
	t183 = t168 * t172 - t169 * t170;
	t145 = t183 * t216;
	t144 = -t151 * qJD(1) + t169 * t195;
	t143 = -t150 * qJD(1) + t168 * t195;
	t132 = t154 * t238;
	t124 = 0.1e1 / t126;
	t120 = 0.1e1 / t122;
	t119 = (-t233 * t163 * t152 + t153 * t188) * t173;
	t116 = t164 * t222 - t221 + (-t152 * t164 + t153 * t215) * t132;
	t115 = t238 * t240 + (qJD(1) * t186 + 0.2e1 * t171 * t181) * t154;
	t1 = [t189 * t216 + (qJD(3) * t186 - t210 * t217) * t154, 0, t115, 0, 0, 0; (t127 * t190 + (t127 * t208 + (-qJD(1) * t119 - t114) * t227) * t124) * t171 + (t128 * t190 * t119 + (((-t123 * t188 - t233 * t208 + t225 * t234) * t152 + (t189 * t218 + t228 + (-t228 + (t234 + t239) * t207) * t154) * t153) * t201 + (t128 * t208 + t163 * t203) * t119 + (t127 + ((-t166 + t167) * t153 * t199 + t233 * t200) * t128) * t163 * qJD(1)) * t124) * t173, 0, 0.2e1 * (-t116 * t227 + t127 * t164) * t173 * t232 + ((t127 * t210 + (qJD(3) * t116 + t114) * t226) * t164 + (t127 * t206 + (t115 * t153 * t171 + t237 * t152 + (qJD(3) * t152 - t123 * t222 + t153 * t209) * t132) * t201 + (-t128 * t210 + t173 * t203) * t116 + ((-t115 + t209) * t152 + ((t132 * t171 - 0.1e1) * qJD(3) + (-t132 + t171) * t123) * t153) * t164 * t226) * t163) * t124, 0, 0, 0; (t134 * t185 - t138 * t224) * t204 + ((t138 * qJD(6) - t143 * t172 + t144 * t170) * t134 + t138 * t118 * t202 + (t185 * t118 + (t185 * qJD(6) + t143 * t170 + t144 * t172) * t184 - t138 * t117) * t135) * t120, 0, (-t134 * t145 - t184 * t223) * t204 + (-t117 * t223 + (-t135 * t145 + t146 * t202) * t118 + (t183 * t134 + t182 * t224) * t164 * t206 + ((t236 * t134 + t235 * t224) * t169 + (-t235 * t134 + t236 * t224) * t168) * t163) * t120, 0, 0, -0.2e1 * t231 - 0.2e1 * (t117 * t135 * t120 - (-t120 * t229 - t135 * t231) * t184) * t184;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end