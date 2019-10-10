% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPP6
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
%   Wie in S6RPRRPP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:18
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:18:23
	% EndTime: 2019-10-10 01:18:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:18:23
	% EndTime: 2019-10-10 01:18:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:18:23
	% EndTime: 2019-10-10 01:18:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:18:23
	% EndTime: 2019-10-10 01:18:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:18:24
	% EndTime: 2019-10-10 01:18:25
	% DurationCPUTime: 0.92s
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
	t191 = t93 * t97;
	t190 = t96 * t93;
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
	t1 = [-0.2e1 * t129 * t120 * t162 + (-qJD(3) * t146 + t120 * t159) * t116, 0, t86, 0, 0, 0; (t169 * t190 + (0.2e1 * t164 + (qJD(1) * t89 + t85) * t191) * t131) * t132 + (t89 * t157 * t131 + (-t89 * t163 + (t89 * t166 + ((t90 * t149 + t194 * t169 + 0.2e1 * t162) * t113 + (t154 * t180 + t131 * t90 + (t124 * t158 + (-t90 + 0.2e1 * t167) * t131) * t116) * t114) * t187) * t131 + (t96 + ((t123 - t126) * t114 * t161 - t132 * t148) * t97) * t171) * t93) * t129, 0, (t170 * t190 + (-0.2e1 * t164 + (-qJD(3) * t88 - t85) * t191) * t129) * t128 + (t88 * t129 * t157 + (t129 * qJD(3) * t96 + (-t114 * t132 * t86 + t152 * t113 + (-qJD(3) * t113 + t114 * t172 + t182 * t90) * t103) * t97 * t176 + (t129 * t166 + t97 * t170) * t88 + ((-t86 - t172) * t113 + (t152 * t103 - t153) * t114) * t128 * t187) * t93) * t131, 0, 0, 0; (-t106 * t111 + t112 * t184) * t165 + (t112 * t155 + t132 * t106 * t144 + t140 * t185 + (t132 * t109 * t145 - t111 * t92 - t112 * t91 - t140 * t183) * t107) * t101, 0, t142 * t165 * t176 + (-t142 * t159 + (t142 * t169 + ((qJD(4) * t106 + t155) * t130 + (-t130 * t91 + (qJD(4) * t109 - t92) * t127) * t107) * t131) * t129) * t101, -0.2e1 * t189 + 0.2e1 * (t91 * t107 * t101 + (-t101 * t188 - t107 * t189) * t109) * t109, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:18:24
	% EndTime: 2019-10-10 01:18:25
	% DurationCPUTime: 0.98s
	% Computational Cost: add. (1161->91), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->95)
	t137 = cos(qJ(1));
	t200 = 0.2e1 * t137;
	t126 = qJ(4) + pkin(9);
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
	% StartTime: 2019-10-10 01:18:24
	% EndTime: 2019-10-10 01:18:25
	% DurationCPUTime: 1.45s
	% Computational Cost: add. (4821->119), mult. (6168->264), div. (1114->15), fcn. (7752->9), ass. (0->112)
	t159 = sin(qJ(1));
	t160 = cos(qJ(3));
	t209 = t159 * t160;
	t152 = qJ(4) + pkin(9);
	t150 = sin(t152);
	t158 = sin(qJ(3));
	t230 = cos(qJ(1));
	t191 = t230 * t158;
	t151 = cos(t152);
	t210 = t159 * t151;
	t139 = t150 * t191 + t210;
	t208 = t160 * t150;
	t127 = atan2(t139, t208);
	t124 = cos(t127);
	t123 = sin(t127);
	t221 = t123 * t139;
	t118 = t124 * t208 + t221;
	t115 = 0.1e1 / t118;
	t138 = t230 * t150 + t158 * t210;
	t133 = 0.1e1 / t138;
	t147 = 0.1e1 / t150;
	t155 = 0.1e1 / t160;
	t116 = 0.1e1 / t118 ^ 2;
	t134 = 0.1e1 / t138 ^ 2;
	t148 = 0.1e1 / t150 ^ 2;
	t211 = t159 * t150;
	t137 = -t230 * t151 + t158 * t211;
	t132 = t137 ^ 2;
	t113 = t116 * t132 + 0.1e1;
	t186 = qJD(1) * t230;
	t204 = qJD(3) * t160;
	t172 = -t158 * t186 - t159 * t204;
	t169 = t230 * qJD(4) - t172;
	t181 = qJD(4) * t158 + qJD(1);
	t121 = t169 * t150 + t181 * t210;
	t224 = t121 * t116;
	t136 = t139 ^ 2;
	t156 = 0.1e1 / t160 ^ 2;
	t214 = t148 * t156;
	t128 = t136 * t214 + 0.1e1;
	t125 = 0.1e1 / t128;
	t202 = qJD(4) * t160;
	t206 = qJD(3) * t158;
	t173 = -t150 * t206 + t151 * t202;
	t193 = t139 * t214;
	t190 = t230 * t160;
	t178 = qJD(3) * t190;
	t179 = t151 * t191;
	t119 = -qJD(4) * t179 - t150 * t178 - t151 * t186 + (qJD(1) * t158 + qJD(4)) * t211;
	t216 = t147 * t155;
	t196 = t119 * t216;
	t107 = (-t173 * t193 - t196) * t125;
	t171 = -t107 * t139 - t173;
	t103 = (-t107 * t208 - t119) * t123 - t171 * t124;
	t117 = t115 * t116;
	t228 = t103 * t117;
	t229 = (-t132 * t228 + t137 * t224) / t113 ^ 2;
	t149 = t147 * t148;
	t154 = t160 ^ 2;
	t157 = t155 / t154;
	t203 = qJD(4) * t151;
	t188 = t156 * t203;
	t227 = (-t119 * t193 + (t148 * t157 * t206 - t149 * t188) * t136) / t128 ^ 2;
	t153 = t159 ^ 2;
	t213 = t153 * t154;
	t195 = t134 * t213;
	t131 = 0.1e1 + t195;
	t170 = -t153 * t158 * t204 + t154 * t159 * t186;
	t122 = t169 * t151 - t181 * t211;
	t223 = t122 * t133 * t134;
	t180 = t213 * t223;
	t226 = (t170 * t134 - t180) / t131 ^ 2;
	t225 = t116 * t137;
	t222 = t123 * t137;
	t220 = t123 * t160;
	t219 = t124 * t137;
	t218 = t124 * t139;
	t217 = t124 * t158;
	t215 = t148 * t151;
	t212 = t158 * t159;
	t207 = qJD(1) * t159;
	t205 = qJD(3) * t159;
	t201 = 0.2e1 * t229;
	t200 = -0.2e1 * t227;
	t199 = 0.2e1 * t226;
	t198 = 0.2e1 * t117 * t137;
	t197 = t116 * t222;
	t194 = t139 * t216;
	t192 = t147 * t156 * t158;
	t189 = t156 * t206;
	t175 = t139 * t192 + t230;
	t114 = t175 * t125;
	t187 = t230 - t114;
	t185 = -0.2e1 * t115 * t229;
	t184 = t116 * t201;
	t183 = 0.2e1 * t155 * t227;
	t182 = -0.2e1 * t137 * t209;
	t177 = t147 * t183;
	t140 = t179 - t211;
	t176 = t139 * t215 - t140 * t147;
	t174 = t134 * t140 * t159 - t230 * t133;
	t168 = t121 * t216 - (t148 * t155 * t203 - t147 * t189) * t137;
	t129 = 0.1e1 / t131;
	t120 = t138 * qJD(1) + t139 * qJD(4) - t151 * t178;
	t111 = 0.1e1 / t113;
	t110 = t176 * t155 * t125;
	t106 = (-t123 + (-t124 * t194 + t123) * t125) * t137;
	t105 = t114 * t218 + (t187 * t220 - t217) * t150;
	t104 = t124 * t151 * t160 + t123 * t140 - (-t123 * t208 + t218) * t110;
	t102 = t175 * t200 + (-t119 * t192 - t207 + (-t148 * t158 * t188 + (0.2e1 * t157 * t158 ^ 2 + t155) * t147 * qJD(3)) * t139) * t125;
	t100 = t176 * t183 + (-t176 * t189 + (t119 * t215 - t120 * t147 + (-t140 * t215 + (0.2e1 * t149 * t151 ^ 2 + t147) * t139) * qJD(4)) * t155) * t125;
	t1 = [-t168 * t125 + t137 * t177, 0, t102, t100, 0, 0; t139 * t185 + (-t119 * t115 + (-t103 * t139 - t106 * t121) * t116) * t111 + (t106 * t184 + (0.2e1 * t106 * t228 + (-t121 * t125 + t121 - (t107 * t125 * t194 + t200) * t137) * t116 * t123 + (-(t139 * t177 - t107) * t225 + (-(t107 + t196) * t137 + t168 * t139) * t116 * t125) * t124) * t111) * t137, 0, t105 * t137 * t184 + (-(t102 * t218 + (-t107 * t221 - t119 * t124) * t114) * t225 + (t103 * t198 - t224) * t105 + (t115 * t209 - (-t114 * t220 + t123 * t190 - t217) * t225) * t203) * t111 + (t185 * t209 + ((-t115 * t205 - (-t187 * qJD(3) + t107) * t197) * t158 + (t115 * t186 + (-t159 * t103 - (-t102 - t207) * t222 - (t187 * t107 - qJD(3)) * t219) * t116) * t160) * t111) * t150, (t104 * t225 - t115 * t138) * t201 + (-t104 * t224 + t122 * t115 + (t104 * t198 - t116 * t138) * t103 - (-t151 * t206 - t150 * t202 + t100 * t139 + t110 * t119 + (t110 * t208 + t140) * t107) * t116 * t219 - (-t120 + (-t100 * t150 - t107 * t151) * t160 - t171 * t110) * t197) * t111, 0, 0; t174 * t160 * t199 + (t174 * t206 + ((-qJD(1) * t133 + 0.2e1 * t140 * t223) * t159 + (t120 * t159 - t230 * t122 - t140 * t186) * t134) * t160) * t129, 0, (t133 * t212 + t151 * t195) * t199 + (0.2e1 * t151 * t180 + t172 * t133 + (qJD(4) * t150 * t213 + t122 * t212 - 0.2e1 * t151 * t170) * t134) * t129, t134 * t182 * t226 + (t182 * t223 + (t121 * t209 + (-t158 * t205 + t160 * t186) * t137) * t134) * t129, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end