% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRP5
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
%   Wie in S6RRRRRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:03
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRP5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:03:55
	% EndTime: 2019-10-10 13:03:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:03:55
	% EndTime: 2019-10-10 13:03:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:03:55
	% EndTime: 2019-10-10 13:03:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:03:56
	% EndTime: 2019-10-10 13:03:57
	% DurationCPUTime: 1.13s
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
	t188 = (t123 * t125 * t155 + (-t123 * t187 - t92 * t151) * t118) / t90 ^ 2;
	t176 = t102 * t106;
	t143 = -qJD(1) * t128 + qJD(3);
	t144 = qJD(3) * t128 - qJD(1);
	t160 = qJD(2) * t129;
	t148 = t125 * t160;
	t87 = -t144 * t166 + (t143 * t126 - t148) * t127;
	t183 = t101 * t102 * t87;
	t167 = t126 * t128;
	t138 = t124 * t167 + t165;
	t86 = t138 * qJD(1) - t107 * qJD(3) + t124 * t148;
	t186 = (-t100 * t183 - t86 * t176) / t99 ^ 2;
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
	t84 = (t189 * t125 * t108 - t109 * t142) * t129;
	t83 = -t192 * t174 + (-t126 + t98) * t128 * t108;
	t81 = t190 * t193 + (qJD(1) * t141 + 0.2e1 * t126 * t137) * t111;
	t1 = [t145 * t171 + (qJD(2) * t141 - t164 * t172) * t111, t81, 0, 0, 0, 0; (-t161 * t184 + (0.2e1 * t157 + (qJD(1) * t84 + t80) * t185) * t125) * t126 + (t84 * t147 * t125 + (-t84 * t155 + (t84 * t159 + ((0.2e1 * t125 * t181 - t85 * t142 - t189 * t161) * t108 + (t118 * t126 * t145 + t125 * t85 + (t117 * t149 - (t85 - 0.2e1 * t162) * t125) * t111) * t109) * t179) * t125 + (-t91 + (-(t119 - t123) * t109 * t153 + t189 * t154) * t92) * t125 * qJD(1)) * t88) * t129, (-t164 * t184 + (-0.2e1 * t157 + (-qJD(2) * t83 - t80) * t185) * t129) * t128 + (t83 * t129 * t147 + (-t91 * t160 - ((-t126 * t81 - t163 * t98) * t109 + (t192 * t85 + t162 - t178) * t108) * t92 * t171 + (t129 * t159 + t92 * t164) * t83 - ((t81 - t163) * t108 + (t85 * t98 + qJD(2) + (-t85 - t178) * t126) * t109) * t128 * t179) * t88) * t125, 0, 0, 0, 0; 0.2e1 * (t101 * t138 + t105 * t176) * t186 + (0.2e1 * t105 * t156 - t144 * t101 * t168 + (t125 * t162 + t140) * t177 + (t138 * t87 + t105 * t86 - t140 * t175 - (qJD(2) * t125 * t127 + t144 * t124) * t106 * t126) * t102) * t96, t128 * t160 * t191 + (-t164 * t191 + (t139 * t158 + ((-qJD(3) * t101 - 0.2e1 * t156) * t127 + (-t127 * t86 + (-qJD(3) * t106 + t87) * t124) * t102) * t96) * t129) * t125, t158 + 0.2e1 * (-t86 * t102 * t96 + (-t102 * t186 - t96 * t183) * t106) * t106, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:03:56
	% EndTime: 2019-10-10 13:03:57
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (1570->96), mult. (2734->204), div. (498->12), fcn. (3199->9), ass. (0->95)
	t152 = sin(qJ(2));
	t145 = t152 ^ 2;
	t154 = cos(qJ(2));
	t148 = 0.1e1 / t154 ^ 2;
	t199 = t145 * t148;
	t153 = sin(qJ(1));
	t217 = 0.2e1 * t153;
	t216 = t152 * t199;
	t151 = qJ(3) + qJ(4);
	t142 = cos(t151);
	t155 = cos(qJ(1));
	t191 = t154 * t155;
	t141 = sin(t151);
	t195 = t153 * t141;
	t131 = t142 * t191 + t195;
	t193 = t153 * t152;
	t135 = atan2(-t193, -t154);
	t134 = cos(t135);
	t133 = sin(t135);
	t180 = t133 * t193;
	t121 = -t134 * t154 - t180;
	t118 = 0.1e1 / t121;
	t125 = 0.1e1 / t131;
	t147 = 0.1e1 / t154;
	t119 = 0.1e1 / t121 ^ 2;
	t126 = 0.1e1 / t131 ^ 2;
	t215 = -0.2e1 * t152;
	t146 = t153 ^ 2;
	t139 = t146 * t199 + 0.1e1;
	t137 = 0.1e1 / t139;
	t214 = t137 - 0.1e1;
	t143 = qJD(3) + qJD(4);
	t192 = t153 * t154;
	t165 = t141 * t192 + t142 * t155;
	t186 = qJD(2) * t155;
	t176 = t152 * t186;
	t109 = t165 * qJD(1) - t131 * t143 + t141 * t176;
	t194 = t153 * t142;
	t130 = t141 * t191 - t194;
	t124 = t130 ^ 2;
	t117 = t124 * t126 + 0.1e1;
	t204 = t126 * t130;
	t170 = -qJD(1) * t154 + t143;
	t171 = t143 * t154 - qJD(1);
	t201 = t141 * t155;
	t110 = -t171 * t201 + (t170 * t153 - t176) * t142;
	t211 = t110 * t125 * t126;
	t213 = (-t109 * t204 - t124 * t211) / t117 ^ 2;
	t189 = qJD(1) * t155;
	t177 = t152 * t189;
	t187 = qJD(2) * t154;
	t188 = qJD(2) * t153;
	t111 = (-(-t153 * t187 - t177) * t147 + t188 * t199) * t137;
	t202 = t134 * t152;
	t105 = (-t111 * t153 + qJD(2)) * t202 + (-t177 + (t111 - t188) * t154) * t133;
	t212 = t105 * t118 * t119;
	t210 = t111 * t133;
	t209 = t111 * t152;
	t208 = t119 * t152;
	t197 = t147 * t152;
	t164 = qJD(2) * (t147 * t216 + t197);
	t168 = t145 * t153 * t189;
	t207 = (t146 * t164 + t148 * t168) / t139 ^ 2;
	t175 = 0.1e1 + t199;
	t123 = t175 * t153 * t137;
	t206 = t123 * t153;
	t205 = t125 * t141;
	t203 = t130 * t142;
	t200 = t145 * t147;
	t150 = t155 ^ 2;
	t198 = t145 * t150;
	t196 = t152 * t155;
	t190 = qJD(1) * t153;
	t114 = t119 * t198 + 0.1e1;
	t185 = 0.2e1 * (-t198 * t212 + (t150 * t152 * t187 - t168) * t119) / t114 ^ 2;
	t184 = -0.2e1 * t213;
	t183 = 0.2e1 * t212;
	t182 = t130 * t211;
	t181 = t119 * t196;
	t179 = t137 * t200;
	t174 = t152 * t185;
	t173 = t207 * t215;
	t172 = t207 * t217;
	t169 = t153 * t179;
	t167 = t175 * t155;
	t166 = t126 * t203 - t205;
	t163 = t152 * t188 + t170 * t155;
	t129 = -t142 * t192 + t201;
	t115 = 0.1e1 / t117;
	t112 = 0.1e1 / t114;
	t108 = (t214 * t152 * t133 - t134 * t169) * t155;
	t107 = -t133 * t192 + t202 + (t133 * t154 - t134 * t193) * t123;
	t106 = -t175 * t172 + (qJD(1) * t167 + t164 * t217) * t137;
	t102 = t184 + 0.2e1 * (-t109 * t115 * t126 + (-t115 * t211 - t126 * t213) * t130) * t130;
	t1 = [t147 * t155 * t173 + (qJD(2) * t167 - t190 * t197) * t137, t106, 0, 0, 0, 0; (t118 * t174 + (-t118 * t187 + (qJD(1) * t108 + t105) * t208) * t112) * t153 + (t119 * t174 * t108 + (-((t111 * t169 + t214 * t187 + t173) * t133 + (t172 * t200 - t209 + (t209 + (t215 - t216) * t188) * t137) * t134) * t181 + (-t119 * t187 + t152 * t183) * t108 + (-t118 + ((-t146 + t150) * t134 * t179 + t214 * t180) * t119) * t152 * qJD(1)) * t112) * t155, (t107 * t208 - t118 * t154) * t155 * t185 + ((-t118 * t190 + (-qJD(2) * t107 - t105) * t155 * t119) * t154 + (-t118 * t186 - (-t106 * t134 * t153 + t133 * t188 + t206 * t210 - t210 + (-qJD(2) * t133 - t134 * t189) * t123) * t181 + (t119 * t190 + t155 * t183) * t107 - ((t106 - t189) * t133 + ((0.1e1 - t206) * qJD(2) + (t123 - t153) * t111) * t134) * t119 * t191) * t152) * t112, 0, 0, 0, 0; 0.2e1 * (t125 * t165 + t129 * t204) * t213 + (0.2e1 * t129 * t182 - t171 * t125 * t194 + t163 * t205 + (-t171 * t130 * t195 + t129 * t109 + t110 * t165 - t163 * t203) * t126) * t115, t166 * t184 * t196 + (t166 * t154 * t186 + (-t166 * t190 + ((-t125 * t143 - 0.2e1 * t182) * t142 + (-t109 * t142 + (-t130 * t143 + t110) * t141) * t126) * t155) * t152) * t115, t102, t102, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:03:56
	% EndTime: 2019-10-10 13:03:57
	% DurationCPUTime: 1.10s
	% Computational Cost: add. (2348->96), mult. (2949->204), div. (516->12), fcn. (3430->9), ass. (0->95)
	t157 = sin(qJ(2));
	t151 = t157 ^ 2;
	t159 = cos(qJ(2));
	t154 = 0.1e1 / t159 ^ 2;
	t204 = t151 * t154;
	t158 = sin(qJ(1));
	t222 = 0.2e1 * t158;
	t221 = t157 * t204;
	t149 = qJ(3) + qJ(4) + qJ(5);
	t147 = cos(t149);
	t160 = cos(qJ(1));
	t196 = t159 * t160;
	t146 = sin(t149);
	t200 = t158 * t146;
	t136 = t147 * t196 + t200;
	t198 = t158 * t157;
	t141 = atan2(-t198, -t159);
	t140 = cos(t141);
	t139 = sin(t141);
	t185 = t139 * t198;
	t126 = -t140 * t159 - t185;
	t123 = 0.1e1 / t126;
	t130 = 0.1e1 / t136;
	t153 = 0.1e1 / t159;
	t124 = 0.1e1 / t126 ^ 2;
	t131 = 0.1e1 / t136 ^ 2;
	t220 = -0.2e1 * t157;
	t152 = t158 ^ 2;
	t144 = t152 * t204 + 0.1e1;
	t142 = 0.1e1 / t144;
	t219 = t142 - 0.1e1;
	t148 = qJD(3) + qJD(4) + qJD(5);
	t197 = t158 * t159;
	t170 = t146 * t197 + t147 * t160;
	t191 = qJD(2) * t160;
	t181 = t157 * t191;
	t114 = t170 * qJD(1) - t136 * t148 + t146 * t181;
	t199 = t158 * t147;
	t135 = t146 * t196 - t199;
	t129 = t135 ^ 2;
	t119 = t129 * t131 + 0.1e1;
	t209 = t131 * t135;
	t175 = -qJD(1) * t159 + t148;
	t176 = t148 * t159 - qJD(1);
	t206 = t146 * t160;
	t115 = -t176 * t206 + (t175 * t158 - t181) * t147;
	t216 = t115 * t130 * t131;
	t218 = (-t114 * t209 - t129 * t216) / t119 ^ 2;
	t194 = qJD(1) * t160;
	t182 = t157 * t194;
	t192 = qJD(2) * t159;
	t193 = qJD(2) * t158;
	t116 = (-(-t158 * t192 - t182) * t153 + t193 * t204) * t142;
	t207 = t140 * t157;
	t110 = (-t116 * t158 + qJD(2)) * t207 + (-t182 + (t116 - t193) * t159) * t139;
	t217 = t110 * t123 * t124;
	t215 = t116 * t139;
	t214 = t116 * t157;
	t213 = t124 * t157;
	t202 = t153 * t157;
	t169 = qJD(2) * (t153 * t221 + t202);
	t173 = t151 * t158 * t194;
	t212 = (t152 * t169 + t154 * t173) / t144 ^ 2;
	t180 = 0.1e1 + t204;
	t128 = t180 * t158 * t142;
	t211 = t128 * t158;
	t210 = t130 * t146;
	t208 = t135 * t147;
	t205 = t151 * t153;
	t156 = t160 ^ 2;
	t203 = t151 * t156;
	t201 = t157 * t160;
	t195 = qJD(1) * t158;
	t122 = t124 * t203 + 0.1e1;
	t190 = 0.2e1 * (-t203 * t217 + (t156 * t157 * t192 - t173) * t124) / t122 ^ 2;
	t189 = -0.2e1 * t218;
	t188 = 0.2e1 * t217;
	t187 = t135 * t216;
	t186 = t124 * t201;
	t184 = t142 * t205;
	t179 = t157 * t190;
	t178 = t212 * t220;
	t177 = t212 * t222;
	t174 = t158 * t184;
	t172 = t180 * t160;
	t171 = t131 * t208 - t210;
	t168 = t157 * t193 + t175 * t160;
	t134 = -t147 * t197 + t206;
	t120 = 0.1e1 / t122;
	t117 = 0.1e1 / t119;
	t113 = (t219 * t157 * t139 - t140 * t174) * t160;
	t112 = -t139 * t197 + t207 + (t139 * t159 - t140 * t198) * t128;
	t111 = -t180 * t177 + (qJD(1) * t172 + t169 * t222) * t142;
	t107 = t189 + 0.2e1 * (-t114 * t117 * t131 + (-t117 * t216 - t131 * t218) * t135) * t135;
	t1 = [t153 * t160 * t178 + (qJD(2) * t172 - t195 * t202) * t142, t111, 0, 0, 0, 0; (t123 * t179 + (-t123 * t192 + (qJD(1) * t113 + t110) * t213) * t120) * t158 + (t124 * t179 * t113 + (-((t116 * t174 + t219 * t192 + t178) * t139 + (t177 * t205 - t214 + (t214 + (t220 - t221) * t193) * t142) * t140) * t186 + (-t124 * t192 + t157 * t188) * t113 + (-t123 + ((-t152 + t156) * t140 * t184 + t219 * t185) * t124) * t157 * qJD(1)) * t120) * t160, (t112 * t213 - t123 * t159) * t160 * t190 + ((-t123 * t195 + (-qJD(2) * t112 - t110) * t160 * t124) * t159 + (-t123 * t191 - (-t111 * t140 * t158 + t139 * t193 + t211 * t215 - t215 + (-qJD(2) * t139 - t140 * t194) * t128) * t186 + (t124 * t195 + t160 * t188) * t112 - ((t111 - t194) * t139 + ((0.1e1 - t211) * qJD(2) + (t128 - t158) * t116) * t140) * t124 * t196) * t157) * t120, 0, 0, 0, 0; 0.2e1 * (t130 * t170 + t134 * t209) * t218 + (0.2e1 * t134 * t187 - t176 * t130 * t199 + t168 * t210 + (-t176 * t135 * t200 + t134 * t114 + t115 * t170 - t168 * t208) * t131) * t117, t171 * t189 * t201 + (t171 * t159 * t191 + (-t171 * t195 + ((-t130 * t148 - 0.2e1 * t187) * t147 + (-t114 * t147 + (-t135 * t148 + t115) * t146) * t131) * t160) * t157) * t117, t107, t107, t107, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:03:56
	% EndTime: 2019-10-10 13:03:57
	% DurationCPUTime: 1.11s
	% Computational Cost: add. (2348->96), mult. (2949->204), div. (516->12), fcn. (3430->9), ass. (0->95)
	t161 = sin(qJ(2));
	t155 = t161 ^ 2;
	t163 = cos(qJ(2));
	t158 = 0.1e1 / t163 ^ 2;
	t208 = t155 * t158;
	t162 = sin(qJ(1));
	t226 = 0.2e1 * t162;
	t225 = t161 * t208;
	t153 = qJ(3) + qJ(4) + qJ(5);
	t151 = cos(t153);
	t164 = cos(qJ(1));
	t200 = t163 * t164;
	t150 = sin(t153);
	t204 = t162 * t150;
	t140 = t151 * t200 + t204;
	t202 = t162 * t161;
	t145 = atan2(-t202, -t163);
	t144 = cos(t145);
	t143 = sin(t145);
	t189 = t143 * t202;
	t130 = -t144 * t163 - t189;
	t127 = 0.1e1 / t130;
	t134 = 0.1e1 / t140;
	t157 = 0.1e1 / t163;
	t128 = 0.1e1 / t130 ^ 2;
	t135 = 0.1e1 / t140 ^ 2;
	t224 = -0.2e1 * t161;
	t156 = t162 ^ 2;
	t148 = t156 * t208 + 0.1e1;
	t146 = 0.1e1 / t148;
	t223 = t146 - 0.1e1;
	t152 = qJD(3) + qJD(4) + qJD(5);
	t201 = t162 * t163;
	t174 = t150 * t201 + t151 * t164;
	t195 = qJD(2) * t164;
	t185 = t161 * t195;
	t118 = t174 * qJD(1) - t140 * t152 + t150 * t185;
	t203 = t162 * t151;
	t139 = t150 * t200 - t203;
	t133 = t139 ^ 2;
	t123 = t133 * t135 + 0.1e1;
	t213 = t135 * t139;
	t179 = -qJD(1) * t163 + t152;
	t180 = t152 * t163 - qJD(1);
	t210 = t150 * t164;
	t119 = -t180 * t210 + (t179 * t162 - t185) * t151;
	t220 = t119 * t134 * t135;
	t222 = (-t118 * t213 - t133 * t220) / t123 ^ 2;
	t198 = qJD(1) * t164;
	t186 = t161 * t198;
	t196 = qJD(2) * t163;
	t197 = qJD(2) * t162;
	t120 = (-(-t162 * t196 - t186) * t157 + t197 * t208) * t146;
	t211 = t144 * t161;
	t114 = (-t120 * t162 + qJD(2)) * t211 + (-t186 + (t120 - t197) * t163) * t143;
	t221 = t114 * t127 * t128;
	t219 = t120 * t143;
	t218 = t120 * t161;
	t217 = t128 * t161;
	t206 = t157 * t161;
	t173 = qJD(2) * (t157 * t225 + t206);
	t177 = t155 * t162 * t198;
	t216 = (t156 * t173 + t158 * t177) / t148 ^ 2;
	t184 = 0.1e1 + t208;
	t132 = t184 * t162 * t146;
	t215 = t132 * t162;
	t214 = t134 * t150;
	t212 = t139 * t151;
	t209 = t155 * t157;
	t160 = t164 ^ 2;
	t207 = t155 * t160;
	t205 = t161 * t164;
	t199 = qJD(1) * t162;
	t126 = t128 * t207 + 0.1e1;
	t194 = 0.2e1 * (-t207 * t221 + (t160 * t161 * t196 - t177) * t128) / t126 ^ 2;
	t193 = -0.2e1 * t222;
	t192 = 0.2e1 * t221;
	t191 = t139 * t220;
	t190 = t128 * t205;
	t188 = t146 * t209;
	t183 = t161 * t194;
	t182 = t216 * t224;
	t181 = t216 * t226;
	t178 = t162 * t188;
	t176 = t184 * t164;
	t175 = t135 * t212 - t214;
	t172 = t161 * t197 + t179 * t164;
	t138 = -t151 * t201 + t210;
	t124 = 0.1e1 / t126;
	t121 = 0.1e1 / t123;
	t117 = (t223 * t161 * t143 - t144 * t178) * t164;
	t116 = -t143 * t201 + t211 + (t143 * t163 - t144 * t202) * t132;
	t115 = -t184 * t181 + (qJD(1) * t176 + t173 * t226) * t146;
	t111 = t193 + 0.2e1 * (-t118 * t121 * t135 + (-t121 * t220 - t135 * t222) * t139) * t139;
	t1 = [t157 * t164 * t182 + (qJD(2) * t176 - t199 * t206) * t146, t115, 0, 0, 0, 0; (t127 * t183 + (-t127 * t196 + (qJD(1) * t117 + t114) * t217) * t124) * t162 + (t128 * t183 * t117 + (-((t120 * t178 + t223 * t196 + t182) * t143 + (t181 * t209 - t218 + (t218 + (t224 - t225) * t197) * t146) * t144) * t190 + (-t128 * t196 + t161 * t192) * t117 + (-t127 + ((-t156 + t160) * t144 * t188 + t223 * t189) * t128) * t161 * qJD(1)) * t124) * t164, (t116 * t217 - t127 * t163) * t164 * t194 + ((-t127 * t199 + (-qJD(2) * t116 - t114) * t164 * t128) * t163 + (-t127 * t195 - (-t115 * t144 * t162 + t143 * t197 + t215 * t219 - t219 + (-qJD(2) * t143 - t144 * t198) * t132) * t190 + (t128 * t199 + t164 * t192) * t116 - ((t115 - t198) * t143 + ((0.1e1 - t215) * qJD(2) + (t132 - t162) * t120) * t144) * t128 * t200) * t161) * t124, 0, 0, 0, 0; 0.2e1 * (t134 * t174 + t138 * t213) * t222 + (0.2e1 * t138 * t191 - t180 * t134 * t203 + t172 * t214 + (-t180 * t139 * t204 + t138 * t118 + t119 * t174 - t172 * t212) * t135) * t121, t175 * t193 * t205 + (t175 * t163 * t195 + (-t175 * t199 + ((-t134 * t152 - 0.2e1 * t191) * t151 + (-t118 * t151 + (-t139 * t152 + t119) * t150) * t135) * t164) * t161) * t121, t111, t111, t111, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end