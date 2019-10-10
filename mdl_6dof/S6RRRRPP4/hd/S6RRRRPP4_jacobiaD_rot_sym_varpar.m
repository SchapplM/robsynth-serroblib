% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP4
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
%   Wie in S6RRRRPP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:25
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:46
	% EndTime: 2019-10-10 12:25:47
	% DurationCPUTime: 1.01s
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
	% StartTime: 2019-10-10 12:25:46
	% EndTime: 2019-10-10 12:25:47
	% DurationCPUTime: 1.06s
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
	% StartTime: 2019-10-10 12:25:46
	% EndTime: 2019-10-10 12:25:47
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (2009->96), mult. (2734->204), div. (498->12), fcn. (3199->9), ass. (0->95)
	t159 = sin(qJ(2));
	t153 = t159 ^ 2;
	t161 = cos(qJ(2));
	t156 = 0.1e1 / t161 ^ 2;
	t206 = t153 * t156;
	t160 = sin(qJ(1));
	t224 = 0.2e1 * t160;
	t223 = t159 * t206;
	t150 = qJ(3) + qJ(4) + pkin(10);
	t149 = cos(t150);
	t162 = cos(qJ(1));
	t198 = t161 * t162;
	t148 = sin(t150);
	t202 = t160 * t148;
	t138 = t149 * t198 + t202;
	t200 = t160 * t159;
	t143 = atan2(-t200, -t161);
	t141 = cos(t143);
	t140 = sin(t143);
	t187 = t140 * t200;
	t128 = -t141 * t161 - t187;
	t125 = 0.1e1 / t128;
	t132 = 0.1e1 / t138;
	t155 = 0.1e1 / t161;
	t126 = 0.1e1 / t128 ^ 2;
	t133 = 0.1e1 / t138 ^ 2;
	t222 = -0.2e1 * t159;
	t154 = t160 ^ 2;
	t146 = t154 * t206 + 0.1e1;
	t144 = 0.1e1 / t146;
	t221 = t144 - 0.1e1;
	t151 = qJD(3) + qJD(4);
	t199 = t160 * t161;
	t172 = t148 * t199 + t149 * t162;
	t193 = qJD(2) * t162;
	t183 = t159 * t193;
	t116 = t172 * qJD(1) - t138 * t151 + t148 * t183;
	t201 = t160 * t149;
	t137 = t148 * t198 - t201;
	t131 = t137 ^ 2;
	t121 = t131 * t133 + 0.1e1;
	t211 = t133 * t137;
	t177 = -qJD(1) * t161 + t151;
	t178 = t151 * t161 - qJD(1);
	t208 = t148 * t162;
	t117 = -t178 * t208 + (t177 * t160 - t183) * t149;
	t218 = t117 * t132 * t133;
	t220 = (-t116 * t211 - t131 * t218) / t121 ^ 2;
	t196 = qJD(1) * t162;
	t184 = t159 * t196;
	t194 = qJD(2) * t161;
	t195 = qJD(2) * t160;
	t118 = (-(-t160 * t194 - t184) * t155 + t195 * t206) * t144;
	t209 = t141 * t159;
	t112 = (-t118 * t160 + qJD(2)) * t209 + (-t184 + (t118 - t195) * t161) * t140;
	t219 = t112 * t125 * t126;
	t217 = t118 * t140;
	t216 = t118 * t159;
	t215 = t126 * t159;
	t204 = t155 * t159;
	t171 = qJD(2) * (t155 * t223 + t204);
	t175 = t153 * t160 * t196;
	t214 = (t154 * t171 + t156 * t175) / t146 ^ 2;
	t182 = 0.1e1 + t206;
	t130 = t182 * t160 * t144;
	t213 = t130 * t160;
	t212 = t132 * t148;
	t210 = t137 * t149;
	t207 = t153 * t155;
	t158 = t162 ^ 2;
	t205 = t153 * t158;
	t203 = t159 * t162;
	t197 = qJD(1) * t160;
	t124 = t126 * t205 + 0.1e1;
	t192 = 0.2e1 * (-t205 * t219 + (t158 * t159 * t194 - t175) * t126) / t124 ^ 2;
	t191 = -0.2e1 * t220;
	t190 = 0.2e1 * t219;
	t189 = t137 * t218;
	t188 = t126 * t203;
	t186 = t144 * t207;
	t181 = t159 * t192;
	t180 = t214 * t222;
	t179 = t214 * t224;
	t176 = t160 * t186;
	t174 = t182 * t162;
	t173 = t133 * t210 - t212;
	t170 = t159 * t195 + t177 * t162;
	t136 = -t149 * t199 + t208;
	t122 = 0.1e1 / t124;
	t119 = 0.1e1 / t121;
	t115 = (t221 * t159 * t140 - t141 * t176) * t162;
	t114 = -t140 * t199 + t209 + (t140 * t161 - t141 * t200) * t130;
	t113 = -t182 * t179 + (qJD(1) * t174 + t171 * t224) * t144;
	t109 = t191 + 0.2e1 * (-t116 * t119 * t133 + (-t119 * t218 - t133 * t220) * t137) * t137;
	t1 = [t155 * t162 * t180 + (qJD(2) * t174 - t197 * t204) * t144, t113, 0, 0, 0, 0; (t125 * t181 + (-t125 * t194 + (qJD(1) * t115 + t112) * t215) * t122) * t160 + (t126 * t181 * t115 + (-((t118 * t176 + t221 * t194 + t180) * t140 + (t179 * t207 - t216 + (t216 + (t222 - t223) * t195) * t144) * t141) * t188 + (-t126 * t194 + t159 * t190) * t115 + (-t125 + ((-t154 + t158) * t141 * t186 + t221 * t187) * t126) * t159 * qJD(1)) * t122) * t162, (t114 * t215 - t125 * t161) * t162 * t192 + ((-t125 * t197 + (-qJD(2) * t114 - t112) * t162 * t126) * t161 + (-t125 * t193 - (-t113 * t141 * t160 + t140 * t195 + t213 * t217 - t217 + (-qJD(2) * t140 - t141 * t196) * t130) * t188 + (t126 * t197 + t162 * t190) * t114 - ((t113 - t196) * t140 + ((0.1e1 - t213) * qJD(2) + (t130 - t160) * t118) * t141) * t126 * t198) * t159) * t122, 0, 0, 0, 0; 0.2e1 * (t132 * t172 + t136 * t211) * t220 + (0.2e1 * t136 * t189 - t178 * t132 * t201 + t170 * t212 + (-t178 * t137 * t202 + t136 * t116 + t117 * t172 - t170 * t210) * t133) * t119, t173 * t191 * t203 + (t173 * t161 * t193 + (-t173 * t197 + ((-t132 * t151 - 0.2e1 * t189) * t149 + (-t116 * t149 + (-t137 * t151 + t117) * t148) * t133) * t162) * t159) * t119, t109, t109, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:46
	% EndTime: 2019-10-10 12:25:47
	% DurationCPUTime: 1.65s
	% Computational Cost: add. (11090->123), mult. (8378->264), div. (1558->15), fcn. (10537->9), ass. (0->116)
	t191 = qJ(3) + qJ(4) + pkin(10);
	t190 = cos(t191);
	t274 = 0.2e1 * t190;
	t189 = sin(t191);
	t199 = cos(qJ(2));
	t200 = cos(qJ(1));
	t247 = t200 * t190;
	t230 = t199 * t247;
	t268 = sin(qJ(1));
	t175 = t268 * t189 + t230;
	t169 = 0.1e1 / t175 ^ 2;
	t198 = sin(qJ(2));
	t193 = t198 ^ 2;
	t197 = t200 ^ 2;
	t252 = t193 * t197;
	t235 = t169 * t252;
	t165 = 0.1e1 + t235;
	t223 = qJD(1) * t268;
	t245 = qJD(2) * t199;
	t208 = t193 * t200 * t223 - t197 * t198 * t245;
	t192 = qJD(3) + qJD(4);
	t244 = qJD(2) * t200;
	t225 = t198 * t244;
	t211 = t199 * t223 + t225;
	t229 = t268 * t192;
	t248 = t200 * t189;
	t154 = (-t192 * t199 + qJD(1)) * t248 + (t229 - t211) * t190;
	t168 = 0.1e1 / t175;
	t263 = t154 * t168 * t169;
	t218 = t252 * t263;
	t273 = (-t208 * t169 - t218) / t165 ^ 2;
	t227 = t268 * t199;
	t171 = t189 * t227 + t247;
	t272 = t171 * t192;
	t250 = t198 * t200;
	t153 = t171 * qJD(1) - t192 * t230 + (t225 - t229) * t189;
	t174 = -t268 * t190 + t199 * t248;
	t186 = 0.1e1 / t189;
	t187 = 0.1e1 / t189 ^ 2;
	t194 = 0.1e1 / t198;
	t195 = 0.1e1 / t198 ^ 2;
	t226 = t195 * t245;
	t254 = t190 * t192;
	t256 = t186 * t194;
	t271 = (t187 * t194 * t254 + t186 * t226) * t174 + t153 * t256;
	t251 = t198 * t189;
	t161 = atan2(-t171, t251);
	t158 = cos(t161);
	t157 = sin(t161);
	t262 = t157 * t171;
	t152 = t158 * t251 - t262;
	t149 = 0.1e1 / t152;
	t150 = 0.1e1 / t152 ^ 2;
	t270 = -0.2e1 * t171;
	t269 = 0.2e1 * t174;
	t166 = t171 ^ 2;
	t255 = t187 * t195;
	t162 = t166 * t255 + 0.1e1;
	t159 = 0.1e1 / t162;
	t253 = t190 * t198;
	t212 = t189 * t245 + t192 * t253;
	t233 = t171 * t255;
	t228 = t268 * t198;
	t216 = qJD(2) * t228;
	t246 = qJD(1) * t200;
	t155 = -t189 * t216 - t192 * t248 - t190 * t223 + (t189 * t246 + t190 * t229) * t199;
	t236 = t155 * t256;
	t141 = (t212 * t233 - t236) * t159;
	t209 = -t141 * t171 + t212;
	t136 = (-t141 * t251 - t155) * t157 + t209 * t158;
	t151 = t149 * t150;
	t267 = t136 * t151;
	t188 = t186 * t187;
	t196 = t194 / t193;
	t231 = t195 * t254;
	t266 = (t155 * t233 + (-t187 * t196 * t245 - t188 * t231) * t166) / t162 ^ 2;
	t265 = t150 * t174;
	t264 = t153 * t150;
	t261 = t157 * t174;
	t260 = t157 * t198;
	t259 = t158 * t171;
	t258 = t158 * t174;
	t257 = t158 * t199;
	t249 = t199 * t200;
	t167 = t174 ^ 2;
	t147 = t150 * t167 + 0.1e1;
	t243 = 0.2e1 * (-t167 * t267 - t174 * t264) / t147 ^ 2;
	t242 = -0.2e1 * t266;
	t241 = 0.2e1 * t273;
	t240 = t151 * t269;
	t239 = t194 * t266;
	t238 = t150 * t261;
	t234 = t171 * t256;
	t232 = t186 * t195 * t199;
	t214 = t171 * t232 + t268;
	t148 = t214 * t159;
	t224 = t268 - t148;
	t222 = t149 * t243;
	t221 = t150 * t243;
	t220 = t250 * t269;
	t219 = t186 * t239;
	t173 = t190 * t227 - t248;
	t215 = t171 * t187 * t190 - t173 * t186;
	t213 = t169 * t173 * t200 - t268 * t168;
	t163 = 0.1e1 / t165;
	t156 = t175 * qJD(1) - t190 * t216 - t272;
	t145 = 0.1e1 / t147;
	t144 = t215 * t194 * t159;
	t140 = (-t157 + (t158 * t234 + t157) * t159) * t174;
	t139 = -t148 * t259 + (t224 * t260 + t257) * t189;
	t138 = t158 * t253 - t157 * t173 + (-t157 * t251 - t259) * t144;
	t137 = t169 * t220 * t273 + (t220 * t263 + (t153 * t250 + (t198 * t223 - t199 * t244) * t174) * t169) * t163;
	t135 = t214 * t242 + (t155 * t232 + t246 + (-t187 * t199 * t231 + (-0.2e1 * t196 * t199 ^ 2 - t194) * t186 * qJD(2)) * t171) * t159;
	t133 = -0.2e1 * t215 * t239 + (-t215 * t226 + ((-t156 - t272) * t186 + (t188 * t254 * t270 + (t173 * t192 + t155) * t187) * t190) * t194) * t159;
	t132 = (t138 * t265 - t149 * t175) * t243 + (t138 * t264 + t154 * t149 + (t138 * t240 - t175 * t150) * t136 - (t190 * t245 - t192 * t251 - t133 * t171 - t144 * t155 + (-t144 * t251 - t173) * t141) * t150 * t258 - (-t156 + (-t133 * t189 - t141 * t190) * t198 - t209 * t144) * t238) * t145;
	t1 = [t271 * t159 + t219 * t269, t135, t133, t133, 0, 0; t171 * t222 + (-t155 * t149 + (t136 * t171 + t140 * t153) * t150) * t145 + (t140 * t221 + (0.2e1 * t140 * t267 + (t153 * t159 - t153 - (-t141 * t159 * t234 + t242) * t174) * t150 * t157 + (-(t219 * t270 - t141) * t265 + (-(t141 + t236) * t174 + t271 * t171) * t150 * t159) * t158) * t145) * t174, t139 * t174 * t221 + (-(-t135 * t259 + (t141 * t262 - t155 * t158) * t148) * t265 + (-t149 * t250 - (-t148 * t260 + t157 * t228 + t257) * t265) * t254 + (t136 * t240 + t264) * t139) * t145 + (t222 * t250 + ((-t149 * t244 - (t224 * qJD(2) - t141) * t238) * t199 + (t149 * t223 + (t200 * t136 - (-t135 + t246) * t261 - (t224 * t141 - qJD(2)) * t258) * t150) * t198) * t145) * t189, t132, t132, 0, 0; t213 * t198 * t241 + (-t213 * t245 + ((qJD(1) * t168 + 0.2e1 * t173 * t263) * t200 + (-t268 * t154 - t156 * t200 + t173 * t223) * t169) * t198) * t163, (t168 * t249 + t190 * t235) * t241 + (t218 * t274 + t211 * t168 + (t189 * t192 * t252 + t154 * t249 + t208 * t274) * t169) * t163, t137, t137, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end