% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRP8
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
%   Wie in S5RRRRP8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRRP8_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP8_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:02:54
	% EndTime: 2019-12-31 22:02:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:02:54
	% EndTime: 2019-12-31 22:02:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:02:54
	% EndTime: 2019-12-31 22:02:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:02:54
	% EndTime: 2019-12-31 22:02:55
	% DurationCPUTime: 0.72s
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
	% StartTime: 2019-12-31 22:02:54
	% EndTime: 2019-12-31 22:02:55
	% DurationCPUTime: 0.76s
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
	t1 = [t147 * t155 * t173 + (qJD(2) * t167 - t190 * t197) * t137, t106, 0, 0, 0; (t118 * t174 + (-t118 * t187 + (qJD(1) * t108 + t105) * t208) * t112) * t153 + (t119 * t174 * t108 + (-((t111 * t169 + t214 * t187 + t173) * t133 + (t172 * t200 - t209 + (t209 + (t215 - t216) * t188) * t137) * t134) * t181 + (-t119 * t187 + t152 * t183) * t108 + (-t118 + ((-t146 + t150) * t134 * t179 + t214 * t180) * t119) * t152 * qJD(1)) * t112) * t155, (t107 * t208 - t118 * t154) * t155 * t185 + ((-t118 * t190 + (-qJD(2) * t107 - t105) * t155 * t119) * t154 + (-t118 * t186 - (-t106 * t134 * t153 + t133 * t188 + t206 * t210 - t210 + (-qJD(2) * t133 - t134 * t189) * t123) * t181 + (t119 * t190 + t155 * t183) * t107 - ((t106 - t189) * t133 + ((0.1e1 - t206) * qJD(2) + (t123 - t153) * t111) * t134) * t119 * t191) * t152) * t112, 0, 0, 0; 0.2e1 * (t125 * t165 + t129 * t204) * t213 + (0.2e1 * t129 * t182 - t171 * t125 * t194 + t163 * t205 + (-t171 * t130 * t195 + t129 * t109 + t110 * t165 - t163 * t203) * t126) * t115, t166 * t184 * t196 + (t166 * t154 * t186 + (-t166 * t190 + ((-t125 * t143 - 0.2e1 * t182) * t142 + (-t109 * t142 + (-t130 * t143 + t110) * t141) * t126) * t155) * t152) * t115, t102, t102, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:02:54
	% EndTime: 2019-12-31 22:02:55
	% DurationCPUTime: 0.77s
	% Computational Cost: add. (1570->96), mult. (2734->204), div. (498->12), fcn. (3199->9), ass. (0->95)
	t153 = sin(qJ(2));
	t146 = t153 ^ 2;
	t155 = cos(qJ(2));
	t149 = 0.1e1 / t155 ^ 2;
	t200 = t146 * t149;
	t154 = sin(qJ(1));
	t218 = 0.2e1 * t154;
	t217 = t153 * t200;
	t152 = qJ(3) + qJ(4);
	t143 = cos(t152);
	t156 = cos(qJ(1));
	t192 = t155 * t156;
	t142 = sin(t152);
	t196 = t154 * t142;
	t132 = t143 * t192 + t196;
	t194 = t154 * t153;
	t136 = atan2(-t194, -t155);
	t135 = cos(t136);
	t134 = sin(t136);
	t181 = t134 * t194;
	t122 = -t135 * t155 - t181;
	t119 = 0.1e1 / t122;
	t126 = 0.1e1 / t132;
	t148 = 0.1e1 / t155;
	t120 = 0.1e1 / t122 ^ 2;
	t127 = 0.1e1 / t132 ^ 2;
	t216 = -0.2e1 * t153;
	t147 = t154 ^ 2;
	t140 = t147 * t200 + 0.1e1;
	t138 = 0.1e1 / t140;
	t215 = t138 - 0.1e1;
	t144 = qJD(3) + qJD(4);
	t193 = t154 * t155;
	t166 = t142 * t193 + t143 * t156;
	t187 = qJD(2) * t156;
	t177 = t153 * t187;
	t110 = t166 * qJD(1) - t132 * t144 + t142 * t177;
	t195 = t154 * t143;
	t131 = t142 * t192 - t195;
	t125 = t131 ^ 2;
	t118 = t125 * t127 + 0.1e1;
	t205 = t127 * t131;
	t171 = -qJD(1) * t155 + t144;
	t172 = t144 * t155 - qJD(1);
	t202 = t142 * t156;
	t111 = -t172 * t202 + (t171 * t154 - t177) * t143;
	t212 = t111 * t126 * t127;
	t214 = (-t110 * t205 - t125 * t212) / t118 ^ 2;
	t190 = qJD(1) * t156;
	t178 = t153 * t190;
	t188 = qJD(2) * t155;
	t189 = qJD(2) * t154;
	t112 = (-(-t154 * t188 - t178) * t148 + t189 * t200) * t138;
	t203 = t135 * t153;
	t106 = (-t112 * t154 + qJD(2)) * t203 + (-t178 + (t112 - t189) * t155) * t134;
	t213 = t106 * t119 * t120;
	t211 = t112 * t134;
	t210 = t112 * t153;
	t209 = t120 * t153;
	t198 = t148 * t153;
	t165 = qJD(2) * (t148 * t217 + t198);
	t169 = t146 * t154 * t190;
	t208 = (t147 * t165 + t149 * t169) / t140 ^ 2;
	t176 = 0.1e1 + t200;
	t124 = t176 * t154 * t138;
	t207 = t124 * t154;
	t206 = t126 * t142;
	t204 = t131 * t143;
	t201 = t146 * t148;
	t151 = t156 ^ 2;
	t199 = t146 * t151;
	t197 = t153 * t156;
	t191 = qJD(1) * t154;
	t115 = t120 * t199 + 0.1e1;
	t186 = 0.2e1 * (-t199 * t213 + (t151 * t153 * t188 - t169) * t120) / t115 ^ 2;
	t185 = -0.2e1 * t214;
	t184 = 0.2e1 * t213;
	t183 = t131 * t212;
	t182 = t120 * t197;
	t180 = t138 * t201;
	t175 = t153 * t186;
	t174 = t208 * t216;
	t173 = t208 * t218;
	t170 = t154 * t180;
	t168 = t176 * t156;
	t167 = t127 * t204 - t206;
	t164 = t153 * t189 + t171 * t156;
	t130 = -t143 * t193 + t202;
	t116 = 0.1e1 / t118;
	t113 = 0.1e1 / t115;
	t109 = (t215 * t153 * t134 - t135 * t170) * t156;
	t108 = -t134 * t193 + t203 + (t134 * t155 - t135 * t194) * t124;
	t107 = -t176 * t173 + (qJD(1) * t168 + t165 * t218) * t138;
	t103 = t185 + 0.2e1 * (-t110 * t116 * t127 + (-t116 * t212 - t127 * t214) * t131) * t131;
	t1 = [t148 * t156 * t174 + (qJD(2) * t168 - t191 * t198) * t138, t107, 0, 0, 0; (t119 * t175 + (-t119 * t188 + (qJD(1) * t109 + t106) * t209) * t113) * t154 + (t120 * t175 * t109 + (-((t112 * t170 + t215 * t188 + t174) * t134 + (t173 * t201 - t210 + (t210 + (t216 - t217) * t189) * t138) * t135) * t182 + (-t120 * t188 + t153 * t184) * t109 + (-t119 + ((-t147 + t151) * t135 * t180 + t215 * t181) * t120) * t153 * qJD(1)) * t113) * t156, (t108 * t209 - t119 * t155) * t156 * t186 + ((-t119 * t191 + (-qJD(2) * t108 - t106) * t156 * t120) * t155 + (-t119 * t187 - (-t107 * t135 * t154 + t134 * t189 + t207 * t211 - t211 + (-qJD(2) * t134 - t135 * t190) * t124) * t182 + (t120 * t191 + t156 * t184) * t108 - ((t107 - t190) * t134 + ((0.1e1 - t207) * qJD(2) + (t124 - t154) * t112) * t135) * t120 * t192) * t153) * t113, 0, 0, 0; 0.2e1 * (t126 * t166 + t130 * t205) * t214 + (0.2e1 * t130 * t183 - t172 * t126 * t195 + t164 * t206 + (-t172 * t131 * t196 + t130 * t110 + t111 * t166 - t164 * t204) * t127) * t116, t167 * t185 * t197 + (t167 * t155 * t187 + (-t167 * t191 + ((-t126 * t144 - 0.2e1 * t183) * t143 + (-t110 * t143 + (-t131 * t144 + t111) * t142) * t127) * t156) * t153) * t116, t103, t103, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end