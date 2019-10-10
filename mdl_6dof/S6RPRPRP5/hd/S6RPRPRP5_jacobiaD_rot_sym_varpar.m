% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRP5
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
%   Wie in S6RPRPRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:36
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:59
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
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:59
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (2618->94), mult. (2519->205), div. (480->12), fcn. (2968->9), ass. (0->94)
	t145 = sin(qJ(1));
	t143 = t145 ^ 2;
	t142 = pkin(9) + qJ(3);
	t138 = sin(t142);
	t133 = t138 ^ 2;
	t140 = cos(t142);
	t135 = 0.1e1 / t140 ^ 2;
	t192 = t133 * t135;
	t128 = t143 * t192 + 0.1e1;
	t132 = t138 * t133;
	t134 = 0.1e1 / t140;
	t189 = t134 * t138;
	t155 = qJD(3) * (t132 * t134 * t135 + t189);
	t146 = cos(qJ(1));
	t180 = qJD(1) * t146;
	t190 = t133 * t145;
	t159 = t180 * t190;
	t198 = (t135 * t159 + t143 * t155) / t128 ^ 2;
	t207 = -0.2e1 * t198;
	t165 = 0.1e1 + t192;
	t206 = t145 * t165;
	t141 = pkin(10) + qJ(5);
	t139 = cos(t141);
	t137 = sin(t141);
	t185 = t145 * t137;
	t186 = t140 * t146;
	t122 = t139 * t186 + t185;
	t184 = t145 * t138;
	t125 = atan2(-t184, -t140);
	t124 = cos(t125);
	t123 = sin(t125);
	t172 = t123 * t184;
	t109 = -t124 * t140 - t172;
	t106 = 0.1e1 / t109;
	t116 = 0.1e1 / t122;
	t107 = 0.1e1 / t109 ^ 2;
	t117 = 0.1e1 / t122 ^ 2;
	t126 = 0.1e1 / t128;
	t205 = t126 - 0.1e1;
	t171 = t126 * t133 * t134;
	t160 = t145 * t171;
	t99 = (t205 * t138 * t123 - t124 * t160) * t146;
	t204 = t107 * t99;
	t178 = qJD(3) * t145;
	t168 = t135 * t178;
	t169 = t138 * t180;
	t100 = (-(-t140 * t178 - t169) * t134 + t133 * t168) * t126;
	t193 = t124 * t138;
	t95 = (-t100 * t145 + qJD(3)) * t193 + (-t169 + (t100 - t178) * t140) * t123;
	t203 = t106 * t107 * t95;
	t156 = t139 * t146 + t140 * t185;
	t177 = qJD(3) * t146;
	t167 = t138 * t177;
	t101 = t156 * qJD(1) - qJD(5) * t122 + t137 * t167;
	t183 = t145 * t139;
	t121 = t137 * t186 - t183;
	t115 = t121 ^ 2;
	t114 = t115 * t117 + 0.1e1;
	t196 = t117 * t121;
	t161 = -qJD(1) * t140 + qJD(5);
	t162 = qJD(5) * t140 - qJD(1);
	t188 = t137 * t146;
	t102 = -t162 * t188 + (t161 * t145 - t167) * t139;
	t201 = t102 * t116 * t117;
	t202 = 0.1e1 / t114 ^ 2 * (-t101 * t196 - t115 * t201);
	t200 = t107 * t138;
	t199 = t107 * t146;
	t197 = t116 * t137;
	t195 = t121 * t139;
	t194 = t123 * t140;
	t144 = t146 ^ 2;
	t191 = t133 * t144;
	t187 = t138 * t146;
	t113 = t126 * t206;
	t182 = t113 - t145;
	t181 = qJD(1) * t145;
	t179 = qJD(3) * t140;
	t105 = t107 * t191 + 0.1e1;
	t176 = 0.2e1 / t105 ^ 2 * (-t191 * t203 + (t138 * t144 * t179 - t159) * t107);
	t175 = 0.2e1 * t203;
	t174 = -0.2e1 * t202;
	t173 = t121 * t201;
	t166 = t113 * t145 - 0.1e1;
	t164 = t138 * t176;
	t163 = t134 * t207;
	t158 = t165 * t146;
	t157 = t117 * t195 - t197;
	t154 = t138 * t178 + t161 * t146;
	t120 = -t140 * t183 + t188;
	t111 = 0.1e1 / t114;
	t103 = 0.1e1 / t105;
	t98 = -t145 * t194 + t193 + (-t124 * t184 + t194) * t113;
	t96 = t206 * t207 + (qJD(1) * t158 + 0.2e1 * t145 * t155) * t126;
	t1 = [t163 * t187 + (qJD(3) * t158 - t181 * t189) * t126, 0, t96, 0, 0, 0; (t106 * t164 + (-t106 * t179 + (qJD(1) * t99 + t95) * t200) * t103) * t145 + (t164 * t204 + (-t179 * t204 + (t99 * t175 + ((-t100 * t160 + 0.2e1 * t138 * t198 - t205 * t179) * t123 + (t163 * t190 + t100 * t138 + (t132 * t168 - (t100 - 0.2e1 * t178) * t138) * t126) * t124) * t199) * t138 + (-t106 + ((-t143 + t144) * t124 * t171 + t205 * t172) * t107) * t138 * qJD(1)) * t103) * t146, 0, (-t106 * t140 + t98 * t200) * t146 * t176 + ((-t106 * t181 + (-qJD(3) * t98 - t95) * t199) * t140 + ((-qJD(3) * t106 + t98 * t175) * t146 + (t98 * t181 - ((t96 - t180) * t123 + (-t166 * qJD(3) + t182 * t100) * t124) * t186 + (-(-t113 * t180 - t145 * t96) * t124 - (-t182 * qJD(3) + t166 * t100) * t123) * t187) * t107) * t138) * t103, 0, 0, 0; 0.2e1 * (t116 * t156 + t120 * t196) * t202 + (0.2e1 * t120 * t173 - t162 * t116 * t183 + t154 * t197 + (-t162 * t121 * t185 + t120 * t101 + t102 * t156 - t154 * t195) * t117) * t111, 0, t157 * t174 * t187 + (t157 * t140 * t177 + (-t157 * t181 + ((-qJD(5) * t116 - 0.2e1 * t173) * t139 + (-t101 * t139 + (-qJD(5) * t121 + t102) * t137) * t117) * t146) * t138) * t111, 0, t174 + 0.2e1 * (-t101 * t111 * t117 + (-t111 * t201 - t117 * t202) * t121) * t121, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:36:00
	% DurationCPUTime: 1.53s
	% Computational Cost: add. (6874->125), mult. (6168->273), div. (1114->15), fcn. (7752->9), ass. (0->114)
	t172 = pkin(10) + qJ(5);
	t168 = sin(t172);
	t173 = pkin(9) + qJ(3);
	t171 = cos(t173);
	t170 = cos(t172);
	t175 = cos(qJ(1));
	t226 = t175 * t170;
	t245 = sin(qJ(1));
	t150 = t245 * t168 + t171 * t226;
	t144 = 0.1e1 / t150 ^ 2;
	t169 = sin(t173);
	t164 = t169 ^ 2;
	t174 = t175 ^ 2;
	t230 = t164 * t174;
	t211 = t144 * t230;
	t140 = 0.1e1 + t211;
	t199 = qJD(1) * t245;
	t223 = qJD(3) * t175;
	t203 = t169 * t223;
	t185 = t171 * t199 + t203;
	t198 = t245 * qJD(5);
	t227 = t175 * t168;
	t129 = (-qJD(5) * t171 + qJD(1)) * t227 + (t198 - t185) * t170;
	t143 = 0.1e1 / t150;
	t240 = t129 * t143 * t144;
	t193 = t230 * t240;
	t204 = qJD(3) * t169 * t174;
	t248 = (-t193 + (-t164 * t175 * t199 + t171 * t204) * t144) / t140 ^ 2;
	t228 = t169 * t175;
	t206 = t245 * t171;
	t146 = t168 * t206 + t226;
	t190 = t168 * t198;
	t220 = qJD(5) * t175;
	t201 = t170 * t220;
	t128 = t146 * qJD(1) + t168 * t203 - t171 * t201 - t190;
	t149 = -t245 * t170 + t171 * t227;
	t161 = 0.1e1 / t168;
	t162 = 0.1e1 / t168 ^ 2;
	t165 = 0.1e1 / t169;
	t166 = 0.1e1 / t169 ^ 2;
	t224 = qJD(3) * t171;
	t205 = t166 * t224;
	t221 = qJD(5) * t170;
	t233 = t161 * t165;
	t247 = (t162 * t165 * t221 + t161 * t205) * t149 + t128 * t233;
	t229 = t169 * t168;
	t136 = atan2(-t146, t229);
	t133 = cos(t136);
	t132 = sin(t136);
	t239 = t132 * t146;
	t127 = t133 * t229 - t239;
	t124 = 0.1e1 / t127;
	t125 = 0.1e1 / t127 ^ 2;
	t246 = 0.2e1 * t149;
	t141 = t146 ^ 2;
	t232 = t162 * t166;
	t137 = t141 * t232 + 0.1e1;
	t134 = 0.1e1 / t137;
	t186 = t168 * t224 + t169 * t221;
	t209 = t146 * t232;
	t207 = t245 * t169;
	t191 = qJD(3) * t207;
	t192 = t170 * t199;
	t225 = qJD(1) * t175;
	t130 = t170 * t198 * t171 - t192 + (t225 * t171 - t191 - t220) * t168;
	t212 = t130 * t233;
	t116 = (t186 * t209 - t212) * t134;
	t183 = -t116 * t146 + t186;
	t112 = (-t116 * t229 - t130) * t132 + t183 * t133;
	t126 = t124 * t125;
	t244 = t112 * t126;
	t163 = t161 * t162;
	t167 = t165 / t164;
	t202 = t166 * t221;
	t243 = (t130 * t209 + (-t162 * t167 * t224 - t163 * t202) * t141) / t137 ^ 2;
	t242 = t125 * t149;
	t241 = t128 * t125;
	t238 = t132 * t149;
	t237 = t132 * t169;
	t236 = t133 * t146;
	t235 = t133 * t149;
	t234 = t133 * t171;
	t231 = t162 * t170;
	t222 = qJD(5) * t168;
	t142 = t149 ^ 2;
	t122 = t125 * t142 + 0.1e1;
	t219 = 0.2e1 * (-t142 * t244 - t149 * t241) / t122 ^ 2;
	t218 = -0.2e1 * t243;
	t217 = 0.2e1 * t248;
	t216 = t126 * t246;
	t215 = t165 * t243;
	t214 = t125 * t238;
	t210 = t146 * t233;
	t208 = t161 * t166 * t171;
	t188 = t146 * t208 + t245;
	t123 = t188 * t134;
	t200 = t245 - t123;
	t197 = t124 * t219;
	t196 = t125 * t219;
	t195 = t228 * t246;
	t194 = t161 * t215;
	t148 = t170 * t206 - t227;
	t189 = t146 * t231 - t148 * t161;
	t187 = t144 * t148 * t175 - t245 * t143;
	t138 = 0.1e1 / t140;
	t131 = t150 * qJD(1) - t170 * t191 - t171 * t190 - t201;
	t120 = 0.1e1 / t122;
	t119 = t189 * t165 * t134;
	t115 = (-t132 + (t133 * t210 + t132) * t134) * t149;
	t114 = -t123 * t236 + (t200 * t237 + t234) * t168;
	t113 = t133 * t169 * t170 - t132 * t148 + (-t132 * t229 - t236) * t119;
	t111 = t188 * t218 + (t130 * t208 + t225 + (-t162 * t171 * t202 + (-0.2e1 * t167 * t171 ^ 2 - t165) * t161 * qJD(3)) * t146) * t134;
	t109 = -0.2e1 * t189 * t215 + (-t189 * t205 + (t130 * t231 - t131 * t161 + (t148 * t231 + (-0.2e1 * t163 * t170 ^ 2 - t161) * t146) * qJD(5)) * t165) * t134;
	t1 = [t247 * t134 + t194 * t246, 0, t111, 0, t109, 0; t146 * t197 + (-t130 * t124 + (t112 * t146 + t115 * t128) * t125) * t120 + (t115 * t196 + (0.2e1 * t115 * t244 + (t128 * t134 - t128 - (-t116 * t134 * t210 + t218) * t149) * t125 * t132 + (-(-0.2e1 * t146 * t194 - t116) * t242 + (-(t116 + t212) * t149 + t247 * t146) * t125 * t134) * t133) * t120) * t149, 0, t114 * t149 * t196 + (-(-t111 * t236 + (t116 * t239 - t130 * t133) * t123) * t242 + (t112 * t216 + t241) * t114 + (-t124 * t228 - (-t123 * t237 + t132 * t207 + t234) * t242) * t221) * t120 + (t197 * t228 + ((-t124 * t223 - (t200 * qJD(3) - t116) * t214) * t171 + (t124 * t199 + (t175 * t112 - (-t111 + t225) * t238 - (t200 * t116 - qJD(3)) * t235) * t125) * t169) * t120) * t168, 0, (t113 * t242 - t124 * t150) * t219 + (t113 * t241 + t129 * t124 + (t113 * t216 - t125 * t150) * t112 - (t170 * t224 - t169 * t222 - t109 * t146 - t119 * t130 + (-t119 * t229 - t148) * t116) * t125 * t235 - (-t131 + (-t109 * t168 - t116 * t170) * t169 - t183 * t119) * t214) * t120, 0; t187 * t169 * t217 + (-t187 * t224 + ((qJD(1) * t143 + 0.2e1 * t148 * t240) * t175 + (-t245 * t129 - t131 * t175 + t148 * t199) * t144) * t169) * t138, 0, (t143 * t171 * t175 + t170 * t211) * t217 + (0.2e1 * t170 * t193 + t185 * t143 + ((t129 * t175 - 0.2e1 * t170 * t204) * t171 + (t174 * t222 + 0.2e1 * t175 * t192) * t164) * t144) * t138, 0, t144 * t195 * t248 + (t195 * t240 + (t128 * t228 + (t169 * t199 - t171 * t223) * t149) * t144) * t138, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end