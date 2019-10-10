% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR8
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
%   Wie in S6RPRRPR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR8_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR8_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:38
	% EndTime: 2019-10-10 01:35:39
	% DurationCPUTime: 0.94s
	% Computational Cost: add. (813->90), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->96)
	t132 = cos(qJ(1));
	t195 = 0.2e1 * t132;
	t167 = qJD(3) * t132;
	t126 = t132 ^ 2;
	t128 = sin(qJ(3));
	t121 = 0.1e1 / t128 ^ 2;
	t131 = cos(qJ(3));
	t125 = t131 ^ 2;
	t179 = t121 * t125;
	t118 = t126 * t179 + 0.1e1;
	t116 = 0.1e1 / t118;
	t120 = 0.1e1 / t128;
	t158 = t121 * t167;
	t129 = sin(qJ(1));
	t171 = qJD(1) * t131;
	t160 = t129 * t171;
	t90 = ((t128 * t167 + t160) * t120 + t125 * t158) * t116;
	t152 = -t90 + t167;
	t173 = t132 * t131;
	t115 = atan2(-t173, t128);
	t113 = sin(t115);
	t114 = cos(t115);
	t99 = -t113 * t173 + t114 * t128;
	t96 = 0.1e1 / t99;
	t127 = sin(qJ(4));
	t175 = t132 * t127;
	t130 = cos(qJ(4));
	t177 = t129 * t130;
	t110 = t128 * t177 + t175;
	t106 = 0.1e1 / t110;
	t107 = 0.1e1 / t110 ^ 2;
	t97 = 0.1e1 / t99 ^ 2;
	t194 = t116 - 0.1e1;
	t123 = t129 ^ 2;
	t170 = qJD(1) * t132;
	t147 = t125 * t129 * t170;
	t169 = qJD(3) * t128;
	t163 = t97 * t169;
	t153 = -t132 * t90 + qJD(3);
	t181 = t114 * t131;
	t85 = t153 * t181 + (t152 * t128 + t160) * t113;
	t192 = t85 * t96 * t97;
	t95 = t123 * t125 * t97 + 0.1e1;
	t193 = (t97 * t147 + (-t125 * t192 - t131 * t163) * t123) / t95 ^ 2;
	t93 = 0.1e1 / t95;
	t191 = t93 * t96;
	t190 = t93 * t97;
	t174 = t132 * t130;
	t178 = t129 * t127;
	t109 = t128 * t178 - t174;
	t105 = t109 ^ 2;
	t104 = t105 * t107 + 0.1e1;
	t184 = t107 * t109;
	t150 = qJD(1) * t128 + qJD(4);
	t143 = t150 * t132;
	t151 = qJD(4) * t128 + qJD(1);
	t145 = t151 * t127;
	t168 = qJD(3) * t131;
	t92 = t130 * t143 + (t130 * t168 - t145) * t129;
	t188 = t106 * t107 * t92;
	t144 = t151 * t130;
	t91 = t129 * t144 + (t129 * t168 + t143) * t127;
	t189 = 0.1e1 / t104 ^ 2 * (-t105 * t188 + t91 * t184);
	t187 = t129 * t97;
	t124 = t131 * t125;
	t141 = qJD(3) * (-t121 * t124 - t131) * t120;
	t186 = (-t121 * t147 + t126 * t141) / t118 ^ 2;
	t185 = t106 * t127;
	t183 = t109 * t130;
	t182 = t113 * t132;
	t180 = t120 * t125;
	t176 = t129 * t131;
	t172 = qJD(1) * t129;
	t166 = -0.2e1 * t192;
	t165 = 0.2e1 * t189;
	t164 = t96 * t193;
	t162 = t131 * t186;
	t161 = t116 * t180;
	t159 = t131 * t170;
	t157 = -0.2e1 * t97 * t193;
	t156 = 0.1e1 + t179;
	t155 = 0.2e1 * t109 * t188;
	t154 = t186 * t195;
	t149 = t132 * t161;
	t148 = t194 * t131 * t113;
	t146 = t156 * t129;
	t142 = t107 * t183 - t185;
	t140 = -t150 * t129 + t131 * t167;
	t112 = t128 * t174 - t178;
	t111 = t128 * t175 + t177;
	t103 = t156 * t132 * t116;
	t101 = 0.1e1 / t104;
	t89 = (-t114 * t149 - t148) * t129;
	t88 = t128 * t182 + t181 + (-t113 * t128 - t114 * t173) * t103;
	t86 = -t156 * t154 + (-qJD(1) * t146 + t141 * t195) * t116;
	t1 = [-0.2e1 * t129 * t120 * t162 + (-qJD(3) * t146 + t120 * t159) * t116, 0, t86, 0, 0, 0; (t169 * t191 + (0.2e1 * t164 + (qJD(1) * t89 + t85) * t190) * t131) * t132 + (t89 * t157 * t131 + (-t89 * t163 + (t89 * t166 + ((t90 * t149 + t194 * t169 + 0.2e1 * t162) * t113 + (t154 * t180 + t131 * t90 + (t124 * t158 + (-t90 + 0.2e1 * t167) * t131) * t116) * t114) * t187) * t131 + (t96 + ((t123 - t126) * t114 * t161 - t132 * t148) * t97) * t171) * t93) * t129, 0, (t170 * t191 + (-0.2e1 * t164 + (-qJD(3) * t88 - t85) * t190) * t129) * t128 + (t88 * t129 * t157 + (t129 * qJD(3) * t96 + (-t114 * t132 * t86 + t152 * t113 + (-qJD(3) * t113 + t114 * t172 + t182 * t90) * t103) * t97 * t176 + (t129 * t166 + t97 * t170) * t88 + ((-t86 - t172) * t113 + (t152 * t103 - t153) * t114) * t128 * t187) * t93) * t131, 0, 0, 0; (-t106 * t111 + t112 * t184) * t165 + (t112 * t155 + t132 * t106 * t144 + t140 * t185 + (t132 * t109 * t145 - t111 * t92 - t112 * t91 - t140 * t183) * t107) * t101, 0, t142 * t165 * t176 + (-t142 * t159 + (t142 * t169 + ((qJD(4) * t106 + t155) * t130 + (-t130 * t91 + (qJD(4) * t109 - t92) * t127) * t107) * t131) * t129) * t101, -0.2e1 * t189 + 0.2e1 * (t91 * t107 * t101 + (-t101 * t188 - t107 * t189) * t109) * t109, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:38
	% EndTime: 2019-10-10 01:35:39
	% DurationCPUTime: 0.97s
	% Computational Cost: add. (1161->91), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->95)
	t137 = cos(qJ(1));
	t200 = 0.2e1 * t137;
	t126 = qJ(4) + pkin(10);
	t124 = sin(t126);
	t134 = sin(qJ(3));
	t125 = cos(t126);
	t135 = sin(qJ(1));
	t181 = t135 * t125;
	t114 = t124 * t137 + t134 * t181;
	t111 = 0.1e1 / t114 ^ 2;
	t179 = t137 * t125;
	t182 = t135 * t124;
	t113 = t134 * t182 - t179;
	t188 = t113 * t125;
	t110 = 0.1e1 / t114;
	t190 = t110 * t124;
	t148 = t111 * t188 - t190;
	t109 = t113 ^ 2;
	t102 = t109 * t111 + 0.1e1;
	t99 = 0.1e1 / t102;
	t199 = t148 * t99;
	t136 = cos(qJ(3));
	t178 = t137 * t136;
	t119 = atan2(-t178, t134);
	t117 = sin(t119);
	t118 = cos(t119);
	t106 = -t117 * t178 + t118 * t134;
	t103 = 0.1e1 / t106;
	t127 = 0.1e1 / t134;
	t104 = 0.1e1 / t106 ^ 2;
	t128 = 0.1e1 / t134 ^ 2;
	t133 = t137 ^ 2;
	t132 = t136 ^ 2;
	t185 = t128 * t132;
	t122 = t133 * t185 + 0.1e1;
	t120 = 0.1e1 / t122;
	t198 = t120 - 0.1e1;
	t130 = t135 ^ 2;
	t184 = t130 * t132;
	t101 = t104 * t184 + 0.1e1;
	t175 = qJD(1) * t137;
	t152 = t132 * t135 * t175;
	t173 = qJD(3) * t136;
	t172 = qJD(3) * t137;
	t162 = t128 * t172;
	t176 = qJD(1) * t136;
	t164 = t135 * t176;
	t96 = ((t134 * t172 + t164) * t127 + t132 * t162) * t120;
	t157 = -t96 + t172;
	t158 = -t137 * t96 + qJD(3);
	t187 = t118 * t136;
	t90 = t158 * t187 + (t157 * t134 + t164) * t117;
	t194 = t103 * t104 * t90;
	t197 = (-t184 * t194 + (-t130 * t134 * t173 + t152) * t104) / t101 ^ 2;
	t189 = t111 * t113;
	t155 = qJD(1) * t134 + qJD(4);
	t146 = t135 * t173 + t155 * t137;
	t156 = qJD(4) * t134 + qJD(1);
	t150 = t124 * t156;
	t95 = t146 * t125 - t135 * t150;
	t193 = t110 * t111 * t95;
	t149 = t125 * t156;
	t94 = t146 * t124 + t135 * t149;
	t196 = 0.1e1 / t102 ^ 2 * (-t109 * t193 + t94 * t189);
	t97 = 0.1e1 / t101;
	t195 = t104 * t97;
	t131 = t136 * t132;
	t147 = qJD(3) * (-t128 * t131 - t136) * t127;
	t191 = (-t128 * t152 + t133 * t147) / t122 ^ 2;
	t186 = t127 * t132;
	t183 = t134 * t137;
	t177 = qJD(1) * t135;
	t174 = qJD(3) * t134;
	t171 = -0.2e1 * t197;
	t170 = 0.2e1 * t196;
	t169 = -0.2e1 * t194;
	t168 = t103 * t197;
	t167 = t97 * t174;
	t166 = t136 * t191;
	t165 = t120 * t186;
	t163 = t136 * t175;
	t161 = 0.1e1 + t185;
	t160 = 0.2e1 * t113 * t193;
	t159 = t191 * t200;
	t154 = t137 * t165;
	t153 = t198 * t136 * t117;
	t151 = t161 * t135;
	t145 = -t155 * t135 + t136 * t172;
	t116 = t134 * t179 - t182;
	t115 = t124 * t183 + t181;
	t108 = t161 * t137 * t120;
	t93 = (-t118 * t154 - t153) * t135;
	t92 = t117 * t183 + t187 + (-t117 * t134 - t118 * t178) * t108;
	t91 = -t161 * t159 + (-qJD(1) * t151 + t147 * t200) * t120;
	t1 = [-0.2e1 * t135 * t127 * t166 + (-qJD(3) * t151 + t127 * t163) * t120, 0, t91, 0, 0, 0; (t103 * t167 + (0.2e1 * t168 + (qJD(1) * t93 + t90) * t195) * t136) * t137 + (t93 * t136 * t97 * t169 + (-t93 * t167 + (t93 * t171 + ((t96 * t154 + t198 * t174 + 0.2e1 * t166) * t117 + (t159 * t186 + t96 * t136 + (t131 * t162 + (-t96 + 0.2e1 * t172) * t136) * t120) * t118) * t97 * t135) * t136) * t104 + (t103 + ((t130 - t133) * t118 * t165 - t137 * t153) * t104) * t97 * t176) * t135, 0, (t103 * t97 * t175 + (-0.2e1 * t168 + (-qJD(3) * t92 - t90) * t195) * t135) * t134 + (((qJD(3) * t103 + t92 * t169) * t135 + (t92 * t175 + ((t108 * t177 - t137 * t91) * t118 + ((t108 * t137 - 0.1e1) * t96 + (-t108 + t137) * qJD(3)) * t117) * t135 * t136) * t104) * t97 + (t92 * t171 + ((-t91 - t177) * t117 + (t157 * t108 - t158) * t118) * t97 * t134) * t104 * t135) * t136, 0, 0, 0; (-t110 * t115 + t116 * t189) * t170 + (t116 * t160 + t137 * t110 * t149 + t145 * t190 + (t137 * t113 * t150 - t115 * t95 - t116 * t94 - t145 * t188) * t111) * t99, 0, -t163 * t199 + (t174 * t199 + (t148 * t170 + ((qJD(4) * t110 + t160) * t125 + (-t125 * t94 + (qJD(4) * t113 - t95) * t124) * t111) * t99) * t136) * t135, -0.2e1 * t196 + 0.2e1 * (t111 * t94 * t99 + (-t111 * t196 - t99 * t193) * t113) * t113, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:38
	% EndTime: 2019-10-10 01:35:39
	% DurationCPUTime: 1.02s
	% Computational Cost: add. (1820->94), mult. (2734->205), div. (498->12), fcn. (3199->9), ass. (0->97)
	t155 = sin(qJ(3));
	t149 = 0.1e1 / t155 ^ 2;
	t157 = cos(qJ(3));
	t153 = t157 ^ 2;
	t202 = t149 * t153;
	t158 = cos(qJ(1));
	t220 = 0.2e1 * t158;
	t219 = t157 * t202;
	t195 = t158 * t157;
	t139 = atan2(-t195, t155);
	t137 = sin(t139);
	t138 = cos(t139);
	t126 = -t137 * t195 + t138 * t155;
	t123 = 0.1e1 / t126;
	t146 = qJ(4) + pkin(10) + qJ(6);
	t144 = sin(t146);
	t145 = cos(t146);
	t156 = sin(qJ(1));
	t198 = t156 * t145;
	t134 = t144 * t158 + t155 * t198;
	t130 = 0.1e1 / t134;
	t148 = 0.1e1 / t155;
	t124 = 0.1e1 / t126 ^ 2;
	t131 = 0.1e1 / t134 ^ 2;
	t154 = t158 ^ 2;
	t142 = t154 * t202 + 0.1e1;
	t140 = 0.1e1 / t142;
	t218 = t140 - 0.1e1;
	t151 = t156 ^ 2;
	t201 = t151 * t153;
	t122 = t124 * t201 + 0.1e1;
	t192 = qJD(1) * t158;
	t173 = t153 * t156 * t192;
	t190 = qJD(3) * t157;
	t193 = qJD(1) * t157;
	t182 = t156 * t193;
	t189 = qJD(3) * t158;
	t116 = ((t155 * t189 + t182) * t148 + t189 * t202) * t140;
	t204 = t138 * t157;
	t110 = (-t116 * t158 + qJD(3)) * t204 + (t182 + (-t116 + t189) * t155) * t137;
	t215 = t110 * t123 * t124;
	t217 = (-t201 * t215 + (-t151 * t155 * t190 + t173) * t124) / t122 ^ 2;
	t147 = qJD(4) + qJD(6);
	t176 = qJD(1) * t155 + t147;
	t166 = t156 * t190 + t176 * t158;
	t177 = t147 * t155 + qJD(1);
	t170 = t145 * t177;
	t114 = t166 * t144 + t156 * t170;
	t196 = t158 * t145;
	t199 = t156 * t144;
	t133 = t155 * t199 - t196;
	t129 = t133 ^ 2;
	t119 = t129 * t131 + 0.1e1;
	t206 = t131 * t133;
	t171 = t144 * t177;
	t115 = t166 * t145 - t156 * t171;
	t214 = t115 * t130 * t131;
	t216 = (t114 * t206 - t129 * t214) / t119 ^ 2;
	t213 = t116 * t137;
	t212 = t116 * t157;
	t211 = t124 * t156;
	t210 = t124 * t157;
	t168 = qJD(3) * (-t157 - t219) * t148;
	t209 = (-t149 * t173 + t154 * t168) / t142 ^ 2;
	t180 = 0.1e1 + t202;
	t128 = t180 * t158 * t140;
	t208 = t128 * t158;
	t207 = t130 * t144;
	t205 = t133 * t145;
	t203 = t148 * t153;
	t200 = t155 * t158;
	t197 = t156 * t157;
	t194 = qJD(1) * t156;
	t191 = qJD(3) * t155;
	t188 = 0.2e1 * t216;
	t187 = -0.2e1 * t215;
	t186 = t157 * t217;
	t185 = t124 * t197;
	t184 = t157 * t209;
	t183 = t140 * t203;
	t181 = t157 * t192;
	t179 = 0.2e1 * t133 * t214;
	t178 = t209 * t220;
	t175 = t158 * t183;
	t174 = t218 * t157 * t137;
	t172 = t180 * t156;
	t169 = t131 * t205 - t207;
	t167 = -t176 * t156 + t157 * t189;
	t136 = t155 * t196 - t199;
	t135 = t144 * t200 + t198;
	t120 = 0.1e1 / t122;
	t117 = 0.1e1 / t119;
	t113 = (-t138 * t175 - t174) * t156;
	t112 = t137 * t200 + t204 + (-t137 * t155 - t138 * t195) * t128;
	t111 = -t180 * t178 + (-qJD(1) * t172 + t168 * t220) * t140;
	t107 = -0.2e1 * t216 + 0.2e1 * (t114 * t117 * t131 + (-t117 * t214 - t131 * t216) * t133) * t133;
	t1 = [-0.2e1 * t156 * t148 * t184 + (-qJD(3) * t172 + t148 * t181) * t140, 0, t111, 0, 0, 0; (0.2e1 * t123 * t186 + (t123 * t191 + (qJD(1) * t113 + t110) * t210) * t120) * t158 + (-0.2e1 * t124 * t186 * t113 + (((t116 * t175 + t218 * t191 + 0.2e1 * t184) * t137 + (t178 * t203 + t212 + (-t212 + (0.2e1 * t157 + t219) * t189) * t140) * t138) * t185 + (-t124 * t191 + t157 * t187) * t113 + (t123 + ((t151 - t154) * t138 * t183 - t158 * t174) * t124) * t193) * t120) * t156, 0, 0.2e1 * (-t112 * t210 - t123 * t155) * t156 * t217 + ((t123 * t192 + (-qJD(3) * t112 - t110) * t211) * t155 + (t156 * qJD(3) * t123 + (-t111 * t138 * t158 + t137 * t189 + t208 * t213 - t213 + (-qJD(3) * t137 + t138 * t194) * t128) * t185 + (t124 * t192 + t156 * t187) * t112 + ((-t111 - t194) * t137 + ((-0.1e1 + t208) * qJD(3) + (-t128 + t158) * t116) * t138) * t155 * t211) * t157) * t120, 0, 0, 0; (-t130 * t135 + t136 * t206) * t188 + (t136 * t179 + t158 * t130 * t170 + t167 * t207 + (t158 * t133 * t171 - t136 * t114 - t135 * t115 - t167 * t205) * t131) * t117, 0, t169 * t188 * t197 + (-t169 * t181 + (t169 * t191 + ((t130 * t147 + t179) * t145 + (-t114 * t145 + (t133 * t147 - t115) * t144) * t131) * t157) * t156) * t117, t107, 0, t107;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end