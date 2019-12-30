% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRRPP8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRPP8_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP8_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_jacobiaD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:53:57
	% EndTime: 2019-12-29 19:53:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:53:57
	% EndTime: 2019-12-29 19:53:57
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:53:57
	% EndTime: 2019-12-29 19:53:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:53:57
	% EndTime: 2019-12-29 19:53:58
	% DurationCPUTime: 1.53s
	% Computational Cost: add. (1002->94), mult. (2519->211), div. (480->12), fcn. (2968->9), ass. (0->92)
	t126 = sin(qJ(1));
	t119 = t126 ^ 2;
	t125 = sin(qJ(2));
	t118 = t125 ^ 2;
	t128 = cos(qJ(2));
	t121 = 0.1e1 / t128 ^ 2;
	t173 = t118 * t121;
	t113 = t119 * t173 + 0.1e1;
	t117 = t125 * t118;
	t120 = 0.1e1 / t128;
	t172 = t120 * t125;
	t137 = qJD(2) * (t117 * t120 * t121 + t172);
	t129 = cos(qJ(1));
	t163 = qJD(1) * t129;
	t151 = t126 * t163;
	t181 = 0.1e1 / t113 ^ 2 * (t119 * t137 + t151 * t173);
	t193 = -0.2e1 * t181;
	t111 = 0.1e1 / t113;
	t146 = 0.1e1 + t173;
	t190 = t126 * t146;
	t98 = t111 * t190;
	t192 = t126 * t98 - 0.1e1;
	t124 = sin(qJ(3));
	t127 = cos(qJ(3));
	t165 = t129 * t127;
	t107 = t126 * t124 + t128 * t165;
	t102 = 0.1e1 / t107 ^ 2;
	t166 = t129 * t124;
	t168 = t126 * t127;
	t106 = t128 * t166 - t168;
	t175 = t106 * t127;
	t101 = 0.1e1 / t107;
	t177 = t101 * t124;
	t139 = t102 * t175 - t177;
	t100 = t106 ^ 2;
	t99 = t100 * t102 + 0.1e1;
	t96 = 0.1e1 / t99;
	t191 = t139 * t96;
	t169 = t126 * t125;
	t110 = atan2(-t169, -t128);
	t109 = cos(t110);
	t108 = sin(t110);
	t154 = t108 * t169;
	t94 = -t109 * t128 - t154;
	t91 = 0.1e1 / t94;
	t92 = 0.1e1 / t94 ^ 2;
	t189 = t111 - 0.1e1;
	t123 = t129 ^ 2;
	t161 = qJD(2) * t128;
	t155 = t92 * t161;
	t150 = t125 * t163;
	t162 = qJD(2) * t126;
	t174 = t109 * t125;
	t149 = t121 * t162;
	t85 = (-(-t126 * t161 - t150) * t120 + t118 * t149) * t111;
	t80 = (-t126 * t85 + qJD(2)) * t174 + (-t150 + (t85 - t162) * t128) * t108;
	t187 = t80 * t91 * t92;
	t90 = t123 * t118 * t92 + 0.1e1;
	t188 = (t123 * t125 * t155 + (-t123 * t187 - t151 * t92) * t118) / t90 ^ 2;
	t176 = t102 * t106;
	t143 = -qJD(1) * t128 + qJD(3);
	t144 = qJD(3) * t128 - qJD(1);
	t160 = qJD(2) * t129;
	t148 = t125 * t160;
	t87 = -t144 * t166 + (t126 * t143 - t148) * t127;
	t183 = t101 * t102 * t87;
	t167 = t126 * t128;
	t138 = t124 * t167 + t165;
	t86 = t138 * qJD(1) - t107 * qJD(3) + t124 * t148;
	t186 = (-t100 * t183 - t176 * t86) / t99 ^ 2;
	t88 = 0.1e1 / t90;
	t185 = t88 * t92;
	t184 = t91 * t88;
	t179 = t129 * t92;
	t178 = qJD(2) * t98;
	t171 = t125 * t129;
	t164 = qJD(1) * t126;
	t159 = 0.2e1 * t187;
	t158 = -0.2e1 * t186;
	t157 = t91 * t188;
	t156 = t106 * t183;
	t153 = t111 * t118 * t120;
	t147 = 0.2e1 * t92 * t188;
	t145 = t120 * t193;
	t142 = t126 * t153;
	t141 = t146 * t129;
	t140 = t143 * t129;
	t105 = -t127 * t167 + t166;
	t84 = (t108 * t125 * t189 - t109 * t142) * t129;
	t83 = -t192 * t174 + (-t126 + t98) * t128 * t108;
	t81 = t190 * t193 + (qJD(1) * t141 + 0.2e1 * t126 * t137) * t111;
	t1 = [t145 * t171 + (qJD(2) * t141 - t164 * t172) * t111, t81, 0, 0, 0; (-t161 * t184 + (0.2e1 * t157 + (qJD(1) * t84 + t80) * t185) * t125) * t126 + (t84 * t147 * t125 + (-t84 * t155 + (t84 * t159 + ((0.2e1 * t125 * t181 - t142 * t85 - t161 * t189) * t108 + (t118 * t126 * t145 + t125 * t85 + (t117 * t149 - (t85 - 0.2e1 * t162) * t125) * t111) * t109) * t179) * t125 + (-t91 + (-(t119 - t123) * t109 * t153 + t189 * t154) * t92) * t125 * qJD(1)) * t88) * t129, (-t164 * t184 + (-0.2e1 * t157 + (-qJD(2) * t83 - t80) * t185) * t129) * t128 + (t83 * t129 * t147 + (-t91 * t160 - ((-t126 * t81 - t163 * t98) * t109 + (t192 * t85 + t162 - t178) * t108) * t92 * t171 + (t129 * t159 + t164 * t92) * t83 - ((t81 - t163) * t108 + (t85 * t98 + qJD(2) + (-t85 - t178) * t126) * t109) * t128 * t179) * t88) * t125, 0, 0, 0; 0.2e1 * (t101 * t138 + t105 * t176) * t186 + (0.2e1 * t105 * t156 - t144 * t101 * t168 + (t125 * t162 + t140) * t177 + (t138 * t87 + t105 * t86 - t140 * t175 - (qJD(2) * t125 * t127 + t124 * t144) * t106 * t126) * t102) * t96, t128 * t160 * t191 + (-t164 * t191 + (t139 * t158 + ((-qJD(3) * t101 - 0.2e1 * t156) * t127 + (-t127 * t86 + (-qJD(3) * t106 + t87) * t124) * t102) * t96) * t129) * t125, t158 + 0.2e1 * (-t86 * t102 * t96 + (-t102 * t186 - t183 * t96) * t106) * t106, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:53:57
	% EndTime: 2019-12-29 19:53:59
	% DurationCPUTime: 2.06s
	% Computational Cost: add. (1786->120), mult. (5988->261), div. (1156->14), fcn. (7606->9), ass. (0->108)
	t154 = sin(qJ(2));
	t158 = cos(qJ(1));
	t229 = t154 * t158;
	t155 = sin(qJ(1));
	t157 = cos(qJ(2));
	t209 = t155 * t157;
	t146 = 0.1e1 / t154;
	t147 = 0.1e1 / t154 ^ 2;
	t148 = t146 * t147;
	t228 = qJD(2) * (0.2e1 * t148 * t157 ^ 2 + t146);
	t153 = sin(qJ(3));
	t156 = cos(qJ(3));
	t208 = t156 * t158;
	t128 = t153 * t209 + t208;
	t199 = qJD(3) * t158;
	t179 = t156 * t199;
	t201 = qJD(3) * t153;
	t180 = t155 * t201;
	t202 = qJD(2) * t158;
	t182 = t154 * t202;
	t112 = t128 * qJD(1) + t153 * t182 - t157 * t179 - t180;
	t207 = t158 * t153;
	t131 = -t155 * t156 + t157 * t207;
	t143 = 0.1e1 / t153;
	t144 = 0.1e1 / t153 ^ 2;
	t203 = qJD(2) * t157;
	t185 = t147 * t203;
	t200 = qJD(3) * t156;
	t217 = t143 * t146;
	t227 = (t144 * t146 * t200 + t143 * t185) * t131 + t112 * t217;
	t211 = t154 * t153;
	t122 = atan2(-t128, t211);
	t117 = cos(t122);
	t116 = sin(t122);
	t223 = t116 * t128;
	t111 = t117 * t211 - t223;
	t108 = 0.1e1 / t111;
	t150 = 0.1e1 / t158;
	t109 = 0.1e1 / t111 ^ 2;
	t151 = 0.1e1 / t158 ^ 2;
	t125 = t128 ^ 2;
	t216 = t144 * t147;
	t123 = t125 * t216 + 0.1e1;
	t118 = 0.1e1 / t123;
	t169 = t153 * t203 + t154 * t200;
	t189 = t128 * t216;
	t210 = t154 * t155;
	t183 = qJD(2) * t210;
	t204 = qJD(1) * t158;
	t205 = qJD(1) * t155;
	t114 = t200 * t209 - t156 * t205 + (t204 * t157 - t183 - t199) * t153;
	t191 = t114 * t217;
	t100 = (t169 * t189 - t191) * t118;
	t167 = -t100 * t128 + t169;
	t96 = (-t100 * t211 - t114) * t116 + t167 * t117;
	t226 = t108 * t109 * t96;
	t145 = t143 * t144;
	t181 = t147 * t200;
	t184 = t148 * t203;
	t225 = (t114 * t189 + (-t144 * t184 - t145 * t181) * t125) / t123 ^ 2;
	t224 = t109 * t131;
	t222 = t116 * t131;
	t221 = t116 * t154;
	t220 = t117 * t128;
	t219 = t117 * t131;
	t218 = t117 * t157;
	t215 = t144 * t156;
	t214 = t147 * t151;
	t213 = t147 * t157;
	t212 = t151 * t155;
	t188 = t143 * t213;
	t172 = t128 * t188 + t155;
	t107 = t172 * t118;
	t206 = -t107 + t155;
	t126 = t131 ^ 2;
	t106 = t109 * t126 + 0.1e1;
	t198 = 0.2e1 / t106 ^ 2 * (-t112 * t224 - t126 * t226);
	t197 = 0.2e1 * t226;
	t196 = -0.2e1 * t225;
	t113 = (-qJD(3) * t157 + qJD(1)) * t207 + (-t182 + (-qJD(1) * t157 + qJD(3)) * t155) * t156;
	t132 = t155 * t153 + t157 * t208;
	t127 = t132 ^ 2;
	t124 = t127 * t214 + 0.1e1;
	t152 = t150 * t151;
	t195 = 0.2e1 * (t113 * t132 * t214 + (t147 * t152 * t205 - t151 * t184) * t127) / t124 ^ 2;
	t194 = t146 * t225;
	t193 = t109 * t222;
	t190 = t128 * t217;
	t187 = t150 * t213;
	t186 = t151 * t205;
	t178 = t108 * t198;
	t177 = t109 * t198;
	t176 = t131 * t197;
	t175 = t146 * t195;
	t173 = t143 * t194;
	t130 = t156 * t209 - t207;
	t171 = t128 * t215 - t130 * t143;
	t170 = t130 * t150 - t132 * t212;
	t120 = 0.1e1 / t124;
	t115 = qJD(1) * t132 - t156 * t183 - t157 * t180 - t179;
	t104 = 0.1e1 / t106;
	t103 = t171 * t146 * t118;
	t99 = (-t116 + (t117 * t190 + t116) * t118) * t131;
	t98 = -t107 * t220 + (t206 * t221 + t218) * t153;
	t97 = t117 * t154 * t156 - t116 * t130 + (-t116 * t211 - t220) * t103;
	t95 = t172 * t196 + (t114 * t188 + t204 + (-t144 * t157 * t181 - t143 * t228) * t128) * t118;
	t93 = -0.2e1 * t171 * t194 + (-t171 * t185 + (t114 * t215 - t115 * t143 + (t130 * t215 + (-0.2e1 * t145 * t156 ^ 2 - t143) * t128) * qJD(3)) * t146) * t118;
	t1 = [t227 * t118 + 0.2e1 * t131 * t173, t95, t93, 0, 0; t128 * t178 + (-t114 * t108 + (t112 * t99 + t128 * t96) * t109) * t104 + (t99 * t177 + (t99 * t197 + (t112 * t118 - t112 - (-t100 * t118 * t190 + t196) * t131) * t109 * t116 + (-(-0.2e1 * t128 * t173 - t100) * t224 + (-(t100 + t191) * t131 + t227 * t128) * t109 * t118) * t117) * t104) * t131, t98 * t131 * t177 + (t98 * t176 + (-(-t95 * t220 + (t100 * t223 - t114 * t117) * t107) * t131 + t98 * t112) * t109 + (-t108 * t229 - (-t107 * t221 + t116 * t210 + t218) * t224) * t200) * t104 + (t178 * t229 + ((-t108 * t202 - (t206 * qJD(2) - t100) * t193) * t157 + (t108 * t205 + (t158 * t96 - (-t95 + t204) * t222 - (t206 * t100 - qJD(2)) * t219) * t109) * t154) * t104) * t153, (-t108 * t132 + t97 * t224) * t198 + (t97 * t176 + t113 * t108 - (-t115 + (-t100 * t156 - t153 * t93) * t154 - t167 * t103) * t193 + (t97 * t112 - t132 * t96 - (t156 * t203 - t154 * t201 - t103 * t114 - t128 * t93 + (-t103 * t211 - t130) * t100) * t219) * t109) * t104, 0, 0; t170 * t175 + (t170 * t185 + (t113 * t212 - t115 * t150 + (-t130 * t212 + (0.2e1 * t152 * t155 ^ 2 + t150) * t132) * qJD(1)) * t146) * t120, (t132 * t187 + t156) * t195 + (-t113 * t187 + t201 + (t150 * t228 - t186 * t213) * t132) * t120, t131 * t150 * t175 + (t112 * t146 * t150 + (-t146 * t186 + t150 * t185) * t131) * t120, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:53:58
	% EndTime: 2019-12-29 19:54:00
	% DurationCPUTime: 1.95s
	% Computational Cost: add. (1786->118), mult. (5988->259), div. (1156->14), fcn. (7606->9), ass. (0->111)
	t148 = sin(qJ(2));
	t152 = cos(qJ(1));
	t226 = t148 * t152;
	t137 = 0.1e1 / t148;
	t138 = 0.1e1 / t148 ^ 2;
	t139 = t137 * t138;
	t151 = cos(qJ(2));
	t225 = qJD(2) * (0.2e1 * t139 * t151 ^ 2 + t137);
	t147 = sin(qJ(3));
	t149 = sin(qJ(1));
	t150 = cos(qJ(3));
	t194 = qJD(3) * t151;
	t176 = t147 * t194;
	t197 = qJD(2) * t152;
	t177 = t148 * t197;
	t200 = qJD(1) * t151;
	t181 = t149 * t200;
	t195 = qJD(3) * t150;
	t199 = qJD(1) * t152;
	t111 = t152 * t176 - t149 * t195 - t147 * t199 + (t177 + t181) * t150;
	t204 = t150 * t152;
	t130 = t149 * t147 + t151 * t204;
	t140 = 0.1e1 / t150;
	t141 = 0.1e1 / t150 ^ 2;
	t198 = qJD(2) * t151;
	t180 = t138 * t198;
	t196 = qJD(3) * t147;
	t214 = t137 * t140;
	t224 = t130 * (-t137 * t141 * t196 + t140 * t180) + t111 * t214;
	t203 = t152 * t147;
	t205 = t149 * t151;
	t127 = t150 * t205 - t203;
	t207 = t148 * t150;
	t120 = atan2(-t127, t207);
	t115 = cos(t120);
	t114 = sin(t120);
	t220 = t114 * t127;
	t109 = t115 * t207 - t220;
	t106 = 0.1e1 / t109;
	t144 = 0.1e1 / t152;
	t107 = 0.1e1 / t109 ^ 2;
	t145 = 0.1e1 / t152 ^ 2;
	t163 = -t148 * t196 + t150 * t198;
	t123 = t127 ^ 2;
	t213 = t138 * t141;
	t121 = t123 * t213 + 0.1e1;
	t116 = 0.1e1 / t121;
	t184 = t127 * t213;
	t208 = t148 * t149;
	t178 = qJD(2) * t208;
	t113 = qJD(1) * t130 - t149 * t176 - t150 * t178 - t152 * t195;
	t186 = t113 * t214;
	t98 = (t163 * t184 - t186) * t116;
	t161 = -t127 * t98 + t163;
	t167 = -t207 * t98 - t113;
	t94 = t114 * t167 + t115 * t161;
	t223 = t106 * t107 * t94;
	t142 = t140 * t141;
	t179 = t139 * t198;
	t222 = 0.1e1 / t121 ^ 2 * (t113 * t184 + (t138 * t142 * t196 - t141 * t179) * t123);
	t221 = t107 * t130;
	t219 = t114 * t130;
	t218 = t114 * t148;
	t217 = t115 * t127;
	t216 = t115 * t130;
	t215 = t115 * t151;
	t212 = t138 * t145;
	t211 = t138 * t151;
	t210 = t141 * t147;
	t209 = t145 * t149;
	t206 = t149 * t150;
	t183 = t140 * t211;
	t166 = t127 * t183 + t149;
	t105 = t166 * t116;
	t202 = -t105 + t149;
	t201 = qJD(1) * t149;
	t125 = t130 ^ 2;
	t104 = t107 * t125 + 0.1e1;
	t193 = 0.2e1 / t104 ^ 2 * (-t111 * t221 - t125 * t223);
	t192 = 0.2e1 * t223;
	t191 = -0.2e1 * t222;
	t168 = -qJD(3) + t200;
	t169 = -qJD(1) + t194;
	t110 = -t169 * t204 + (t149 * t168 + t177) * t147;
	t129 = -t151 * t203 + t206;
	t124 = t129 ^ 2;
	t122 = t124 * t212 + 0.1e1;
	t146 = t144 * t145;
	t190 = 0.2e1 * (t129 * t110 * t212 + (t138 * t146 * t201 - t145 * t179) * t124) / t122 ^ 2;
	t189 = t137 * t222;
	t188 = t107 * t219;
	t185 = t127 * t214;
	t182 = t144 * t211;
	t175 = t106 * t193;
	t174 = t107 * t193;
	t173 = t130 * t192;
	t172 = t137 * t190;
	t170 = t140 * t189;
	t126 = t147 * t205 + t204;
	t165 = -t126 * t140 + t127 * t210;
	t164 = -t126 * t144 - t129 * t209;
	t118 = 0.1e1 / t122;
	t112 = t169 * t206 + (t152 * t168 - t178) * t147;
	t102 = 0.1e1 / t104;
	t101 = t165 * t137 * t116;
	t97 = (-t114 + (t115 * t185 + t114) * t116) * t130;
	t96 = -t105 * t217 + (t202 * t218 + t215) * t150;
	t95 = -t115 * t148 * t147 + t114 * t126 - (-t114 * t207 - t217) * t101;
	t93 = t166 * t191 + (t113 * t183 + t199 + (-t140 * t225 + t176 * t213) * t127) * t116;
	t91 = 0.2e1 * t165 * t189 + (t165 * t180 + (-t113 * t210 + t112 * t140 + (t126 * t210 + (-0.2e1 * t142 * t147 ^ 2 - t140) * t127) * qJD(3)) * t137) * t116;
	t1 = [t224 * t116 + 0.2e1 * t130 * t170, t93, t91, 0, 0; t127 * t175 + (-t113 * t106 + (t111 * t97 + t127 * t94) * t107) * t102 + (t97 * t174 + (t97 * t192 + (t111 * t116 - t111 - (-t116 * t185 * t98 + t191) * t130) * t107 * t114 + (-(-0.2e1 * t127 * t170 - t98) * t221 + (-(t98 + t186) * t130 + t224 * t127) * t107 * t116) * t115) * t102) * t130, t96 * t130 * t174 + (t96 * t173 + (-(-t217 * t93 + (-t113 * t115 + t220 * t98) * t105) * t130 + t96 * t111) * t107 + (t106 * t226 - (t105 * t218 - t114 * t208 - t215) * t221) * t196) * t102 + (t175 * t226 + ((-t106 * t197 - (qJD(2) * t202 - t98) * t188) * t151 + (t106 * t201 + (t152 * t94 - (-t93 + t199) * t219 - (t202 * t98 - qJD(2)) * t216) * t107) * t148) * t102) * t150, (-t106 * t129 + t221 * t95) * t193 + (t95 * t173 + t110 * t106 - (t112 + (t147 * t98 - t150 * t91) * t148 + t161 * t101) * t188 + (t95 * t111 - t129 * t94 - (-t101 * t167 + t126 * t98 - t127 * t91 - t147 * t198 - t148 * t195) * t216) * t107) * t102, 0, 0; t164 * t172 + (t164 * t180 + (t110 * t209 + t112 * t144 + (t126 * t209 + (0.2e1 * t146 * t149 ^ 2 + t144) * t129) * qJD(1)) * t137) * t118, (t129 * t182 - t147) * t190 + (-t110 * t182 + t195 + (t144 * t225 - t181 * t212) * t129) * t118, t130 * t144 * t172 + (t111 * t137 * t144 + (-t137 * t145 * t201 + t144 * t180) * t130) * t118, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end