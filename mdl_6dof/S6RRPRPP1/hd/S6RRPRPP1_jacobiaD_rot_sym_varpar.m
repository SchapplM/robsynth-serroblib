% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPP1
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
%   Wie in S6RRPRPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:55
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPP1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPP1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:39
	% EndTime: 2019-10-10 09:55:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:39
	% EndTime: 2019-10-10 09:55:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:39
	% EndTime: 2019-10-10 09:55:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:39
	% EndTime: 2019-10-10 09:55:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:40
	% EndTime: 2019-10-10 09:55:41
	% DurationCPUTime: 1.01s
	% Computational Cost: add. (2270->94), mult. (2519->210), div. (480->12), fcn. (2968->9), ass. (0->93)
	t139 = sin(qJ(1));
	t136 = t139 ^ 2;
	t135 = qJ(2) + pkin(9);
	t133 = sin(t135);
	t129 = t133 ^ 2;
	t134 = cos(t135);
	t131 = 0.1e1 / t134 ^ 2;
	t188 = t129 * t131;
	t124 = t136 * t188 + 0.1e1;
	t128 = t133 * t129;
	t130 = 0.1e1 / t134;
	t185 = t130 * t133;
	t149 = qJD(2) * (t128 * t130 * t131 + t185);
	t141 = cos(qJ(1));
	t177 = qJD(1) * t141;
	t186 = t129 * t139;
	t154 = t177 * t186;
	t194 = (t131 * t154 + t136 * t149) / t124 ^ 2;
	t204 = -0.2e1 * t194;
	t161 = 0.1e1 + t188;
	t203 = t139 * t161;
	t140 = cos(qJ(4));
	t179 = t141 * t140;
	t138 = sin(qJ(4));
	t182 = t139 * t138;
	t120 = t134 * t179 + t182;
	t183 = t139 * t133;
	t121 = atan2(-t183, -t134);
	t116 = cos(t121);
	t115 = sin(t121);
	t168 = t115 * t183;
	t105 = -t116 * t134 - t168;
	t102 = 0.1e1 / t105;
	t112 = 0.1e1 / t120;
	t103 = 0.1e1 / t105 ^ 2;
	t113 = 0.1e1 / t120 ^ 2;
	t122 = 0.1e1 / t124;
	t202 = t122 - 0.1e1;
	t137 = t141 ^ 2;
	t176 = qJD(2) * t134;
	t187 = t129 * t137;
	t175 = qJD(2) * t139;
	t163 = t131 * t175;
	t164 = t133 * t177;
	t96 = (-(-t134 * t175 - t164) * t130 + t129 * t163) * t122;
	t158 = t96 - t175;
	t159 = -t139 * t96 + qJD(2);
	t190 = t116 * t133;
	t91 = t159 * t190 + (t158 * t134 - t164) * t115;
	t199 = t102 * t103 * t91;
	t99 = t103 * t187 + 0.1e1;
	t201 = (-t187 * t199 + (t133 * t137 * t176 - t154) * t103) / t99 ^ 2;
	t97 = 0.1e1 / t99;
	t200 = t103 * t97;
	t150 = t134 * t182 + t179;
	t174 = qJD(2) * t141;
	t162 = t133 * t174;
	t100 = t150 * qJD(1) - t120 * qJD(4) + t138 * t162;
	t180 = t141 * t138;
	t181 = t139 * t140;
	t119 = t134 * t180 - t181;
	t111 = t119 ^ 2;
	t110 = t111 * t113 + 0.1e1;
	t192 = t113 * t119;
	t156 = -qJD(1) * t134 + qJD(4);
	t157 = qJD(4) * t134 - qJD(1);
	t101 = -t157 * t180 + (t156 * t139 - t162) * t140;
	t196 = t101 * t112 * t113;
	t198 = 0.1e1 / t110 ^ 2 * (-t100 * t192 - t111 * t196);
	t193 = t112 * t138;
	t191 = t115 * t134;
	t189 = t119 * t140;
	t184 = t133 * t141;
	t178 = qJD(1) * t139;
	t173 = 0.2e1 * t201;
	t172 = 0.2e1 * t199;
	t171 = -0.2e1 * t198;
	t170 = t102 * t201;
	t169 = t97 * t176;
	t167 = t119 * t196;
	t166 = t122 * t129 * t130;
	t160 = t130 * t204;
	t155 = t139 * t166;
	t153 = t161 * t141;
	t152 = t156 * t141;
	t151 = t113 * t189 - t193;
	t118 = -t134 * t181 + t180;
	t108 = 0.1e1 / t110;
	t107 = t122 * t203;
	t95 = (t202 * t133 * t115 - t116 * t155) * t141;
	t93 = -t139 * t191 + t190 + (-t116 * t183 + t191) * t107;
	t92 = t203 * t204 + (qJD(1) * t153 + 0.2e1 * t139 * t149) * t122;
	t1 = [t160 * t184 + (qJD(2) * t153 - t178 * t185) * t122, t92, 0, 0, 0, 0; (-t102 * t169 + (0.2e1 * t170 + (qJD(1) * t95 + t91) * t200) * t133) * t139 + ((-t95 * t169 + (t95 * t173 + ((0.2e1 * t133 * t194 - t96 * t155 - t202 * t176) * t115 + (t160 * t186 + t133 * t96 + (t128 * t163 - (t96 - 0.2e1 * t175) * t133) * t122) * t116) * t97 * t141) * t133) * t103 + (t95 * t172 + (-t102 + ((-t136 + t137) * t116 * t166 + t202 * t168) * t103) * qJD(1)) * t133 * t97) * t141, (-t102 * t97 * t178 + (-0.2e1 * t170 + (-qJD(2) * t93 - t91) * t200) * t141) * t134 + (((-qJD(2) * t102 + t93 * t172) * t141 + (t93 * t178 + (-(-t107 * t177 - t139 * t92) * t116 - ((t107 * t139 - 0.1e1) * t96 + (-t107 + t139) * qJD(2)) * t115) * t184) * t103) * t97 + (t93 * t173 - ((t92 - t177) * t115 + (t158 * t107 + t159) * t116) * t97 * t134) * t103 * t141) * t133, 0, 0, 0, 0; 0.2e1 * (t112 * t150 + t118 * t192) * t198 + (0.2e1 * t118 * t167 - t157 * t112 * t181 + (t133 * t175 + t152) * t193 + (t118 * t100 + t150 * t101 - t152 * t189 - (qJD(2) * t133 * t140 + t157 * t138) * t119 * t139) * t113) * t108, t151 * t171 * t184 + (t151 * t134 * t174 + (-t151 * t178 + ((-qJD(4) * t112 - 0.2e1 * t167) * t140 + (-t100 * t140 + (-qJD(4) * t119 + t101) * t138) * t113) * t141) * t133) * t108, 0, t171 + 0.2e1 * (-t100 * t113 * t108 + (-t108 * t196 - t113 * t198) * t119) * t119, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:40
	% EndTime: 2019-10-10 09:55:41
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (2618->95), mult. (2519->203), div. (480->12), fcn. (2968->9), ass. (0->93)
	t147 = qJ(2) + pkin(9);
	t143 = sin(t147);
	t138 = t143 ^ 2;
	t145 = cos(t147);
	t140 = 0.1e1 / t145 ^ 2;
	t195 = t138 * t140;
	t150 = sin(qJ(1));
	t212 = 0.2e1 * t150;
	t211 = t143 * t195;
	t148 = t150 ^ 2;
	t133 = t148 * t195 + 0.1e1;
	t131 = 0.1e1 / t133;
	t139 = 0.1e1 / t145;
	t151 = cos(qJ(1));
	t185 = qJD(1) * t151;
	t173 = t143 * t185;
	t183 = qJD(2) * t150;
	t105 = (-(-t145 * t183 - t173) * t139 + t183 * t195) * t131;
	t210 = t105 - t183;
	t146 = qJ(4) + pkin(10);
	t144 = cos(t146);
	t142 = sin(t146);
	t189 = t150 * t142;
	t190 = t145 * t151;
	t127 = t144 * t190 + t189;
	t188 = t150 * t143;
	t130 = atan2(-t188, -t145);
	t129 = cos(t130);
	t128 = sin(t130);
	t176 = t128 * t188;
	t114 = -t129 * t145 - t176;
	t111 = 0.1e1 / t114;
	t121 = 0.1e1 / t127;
	t112 = 0.1e1 / t114 ^ 2;
	t122 = 0.1e1 / t127 ^ 2;
	t209 = -0.2e1 * t143;
	t208 = t131 - 0.1e1;
	t197 = t129 * t143;
	t100 = (-t105 * t150 + qJD(2)) * t197 + (t210 * t145 - t173) * t128;
	t207 = t100 * t111 * t112;
	t161 = t144 * t151 + t145 * t189;
	t182 = qJD(2) * t151;
	t172 = t143 * t182;
	t106 = t161 * qJD(1) - t127 * qJD(4) + t142 * t172;
	t187 = t150 * t144;
	t126 = t142 * t190 - t187;
	t120 = t126 ^ 2;
	t119 = t120 * t122 + 0.1e1;
	t200 = t122 * t126;
	t166 = -qJD(1) * t145 + qJD(4);
	t167 = qJD(4) * t145 - qJD(1);
	t192 = t142 * t151;
	t107 = -t167 * t192 + (t166 * t150 - t172) * t144;
	t204 = t107 * t121 * t122;
	t206 = (-t106 * t200 - t120 * t204) / t119 ^ 2;
	t205 = t105 * t143;
	t203 = t112 * t143;
	t193 = t139 * t143;
	t160 = qJD(2) * (t139 * t211 + t193);
	t164 = t138 * t150 * t185;
	t202 = (t140 * t164 + t148 * t160) / t133 ^ 2;
	t201 = t121 * t142;
	t199 = t126 * t144;
	t198 = t128 * t150;
	t196 = t138 * t139;
	t149 = t151 ^ 2;
	t194 = t138 * t149;
	t191 = t143 * t151;
	t186 = qJD(1) * t150;
	t184 = qJD(2) * t145;
	t110 = t112 * t194 + 0.1e1;
	t181 = 0.2e1 / t110 ^ 2 * (-t194 * t207 + (t143 * t149 * t184 - t164) * t112);
	t180 = 0.2e1 * t207;
	t179 = -0.2e1 * t206;
	t178 = t126 * t204;
	t177 = t112 * t191;
	t175 = t131 * t196;
	t171 = 0.1e1 + t195;
	t170 = t143 * t181;
	t169 = t202 * t209;
	t168 = t202 * t212;
	t165 = t150 * t175;
	t163 = t171 * t151;
	t162 = t122 * t199 - t201;
	t159 = t143 * t183 + t166 * t151;
	t125 = -t145 * t187 + t192;
	t118 = t171 * t150 * t131;
	t116 = 0.1e1 / t119;
	t108 = 0.1e1 / t110;
	t104 = (t208 * t143 * t128 - t129 * t165) * t151;
	t103 = -t145 * t198 + t197 + (t128 * t145 - t129 * t188) * t118;
	t101 = -t171 * t168 + (qJD(1) * t163 + t160 * t212) * t131;
	t1 = [t139 * t151 * t169 + (qJD(2) * t163 - t186 * t193) * t131, t101, 0, 0, 0, 0; (t111 * t170 + (-t111 * t184 + (qJD(1) * t104 + t100) * t203) * t108) * t150 + (t112 * t170 * t104 + (-((t105 * t165 + t208 * t184 + t169) * t128 + (t168 * t196 - t205 + (t205 + (t209 - t211) * t183) * t131) * t129) * t177 + (-t112 * t184 + t143 * t180) * t104 + (-t111 + ((-t148 + t149) * t129 * t175 + t208 * t176) * t112) * t143 * qJD(1)) * t108) * t151, (t103 * t203 - t111 * t145) * t151 * t181 + ((-t111 * t186 + (-qJD(2) * t103 - t100) * t151 * t112) * t145 + (-t111 * t182 - (-t101 * t129 * t150 - t210 * t128 + (-qJD(2) * t128 + t105 * t198 - t129 * t185) * t118) * t177 + (t112 * t186 + t151 * t180) * t103 - ((t101 - t185) * t128 + ((-t118 * t150 + 0.1e1) * qJD(2) + (t118 - t150) * t105) * t129) * t112 * t190) * t143) * t108, 0, 0, 0, 0; 0.2e1 * (t121 * t161 + t125 * t200) * t206 + (0.2e1 * t125 * t178 - t167 * t121 * t187 + t159 * t201 + (-t167 * t126 * t189 + t125 * t106 + t107 * t161 - t159 * t199) * t122) * t116, t162 * t179 * t191 + (t162 * t145 * t182 + (-t162 * t186 + ((-qJD(4) * t121 - 0.2e1 * t178) * t144 + (-t106 * t144 + (-qJD(4) * t126 + t107) * t142) * t122) * t151) * t143) * t116, 0, t179 + 0.2e1 * (-t106 * t116 * t122 + (-t116 * t204 - t122 * t206) * t126) * t126, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:40
	% EndTime: 2019-10-10 09:55:41
	% DurationCPUTime: 1.52s
	% Computational Cost: add. (6874->125), mult. (6168->273), div. (1114->15), fcn. (7752->9), ass. (0->114)
	t178 = qJ(4) + pkin(10);
	t174 = sin(t178);
	t179 = qJ(2) + pkin(9);
	t177 = cos(t179);
	t176 = cos(t178);
	t181 = cos(qJ(1));
	t232 = t181 * t176;
	t251 = sin(qJ(1));
	t156 = t251 * t174 + t177 * t232;
	t150 = 0.1e1 / t156 ^ 2;
	t175 = sin(t179);
	t170 = t175 ^ 2;
	t180 = t181 ^ 2;
	t236 = t170 * t180;
	t217 = t150 * t236;
	t146 = 0.1e1 + t217;
	t205 = qJD(1) * t251;
	t229 = qJD(2) * t181;
	t209 = t175 * t229;
	t191 = t177 * t205 + t209;
	t204 = t251 * qJD(4);
	t233 = t181 * t174;
	t135 = (-qJD(4) * t177 + qJD(1)) * t233 + (t204 - t191) * t176;
	t149 = 0.1e1 / t156;
	t246 = t135 * t149 * t150;
	t199 = t236 * t246;
	t210 = qJD(2) * t175 * t180;
	t254 = (-t199 + (-t170 * t181 * t205 + t177 * t210) * t150) / t146 ^ 2;
	t234 = t175 * t181;
	t212 = t251 * t177;
	t152 = t174 * t212 + t232;
	t196 = t174 * t204;
	t226 = qJD(4) * t181;
	t207 = t176 * t226;
	t134 = t152 * qJD(1) + t174 * t209 - t177 * t207 - t196;
	t155 = -t251 * t176 + t177 * t233;
	t167 = 0.1e1 / t174;
	t168 = 0.1e1 / t174 ^ 2;
	t171 = 0.1e1 / t175;
	t172 = 0.1e1 / t175 ^ 2;
	t230 = qJD(2) * t177;
	t211 = t172 * t230;
	t227 = qJD(4) * t176;
	t239 = t167 * t171;
	t253 = (t168 * t171 * t227 + t167 * t211) * t155 + t134 * t239;
	t235 = t175 * t174;
	t142 = atan2(-t152, t235);
	t139 = cos(t142);
	t138 = sin(t142);
	t245 = t138 * t152;
	t133 = t139 * t235 - t245;
	t130 = 0.1e1 / t133;
	t131 = 0.1e1 / t133 ^ 2;
	t252 = 0.2e1 * t155;
	t147 = t152 ^ 2;
	t238 = t168 * t172;
	t143 = t147 * t238 + 0.1e1;
	t140 = 0.1e1 / t143;
	t192 = t174 * t230 + t175 * t227;
	t215 = t152 * t238;
	t213 = t251 * t175;
	t197 = qJD(2) * t213;
	t198 = t176 * t205;
	t231 = qJD(1) * t181;
	t136 = t176 * t204 * t177 - t198 + (t231 * t177 - t197 - t226) * t174;
	t218 = t136 * t239;
	t122 = (t192 * t215 - t218) * t140;
	t189 = -t122 * t152 + t192;
	t118 = (-t122 * t235 - t136) * t138 + t189 * t139;
	t132 = t130 * t131;
	t250 = t118 * t132;
	t169 = t167 * t168;
	t173 = t171 / t170;
	t208 = t172 * t227;
	t249 = (t136 * t215 + (-t168 * t173 * t230 - t169 * t208) * t147) / t143 ^ 2;
	t248 = t131 * t155;
	t247 = t134 * t131;
	t244 = t138 * t155;
	t243 = t138 * t175;
	t242 = t139 * t152;
	t241 = t139 * t155;
	t240 = t139 * t177;
	t237 = t168 * t176;
	t228 = qJD(4) * t174;
	t148 = t155 ^ 2;
	t128 = t131 * t148 + 0.1e1;
	t225 = 0.2e1 * (-t148 * t250 - t155 * t247) / t128 ^ 2;
	t224 = -0.2e1 * t249;
	t223 = 0.2e1 * t254;
	t222 = t132 * t252;
	t221 = t171 * t249;
	t220 = t131 * t244;
	t216 = t152 * t239;
	t214 = t167 * t172 * t177;
	t194 = t152 * t214 + t251;
	t129 = t194 * t140;
	t206 = t251 - t129;
	t203 = t130 * t225;
	t202 = t131 * t225;
	t201 = t234 * t252;
	t200 = t167 * t221;
	t154 = t176 * t212 - t233;
	t195 = t152 * t237 - t154 * t167;
	t193 = t150 * t154 * t181 - t251 * t149;
	t144 = 0.1e1 / t146;
	t137 = t156 * qJD(1) - t176 * t197 - t177 * t196 - t207;
	t126 = 0.1e1 / t128;
	t125 = t195 * t171 * t140;
	t121 = (-t138 + (t139 * t216 + t138) * t140) * t155;
	t120 = -t129 * t242 + (t206 * t243 + t240) * t174;
	t119 = t139 * t175 * t176 - t138 * t154 + (-t138 * t235 - t242) * t125;
	t117 = t194 * t224 + (t136 * t214 + t231 + (-t168 * t177 * t208 + (-0.2e1 * t173 * t177 ^ 2 - t171) * t167 * qJD(2)) * t152) * t140;
	t115 = -0.2e1 * t195 * t221 + (-t195 * t211 + (t136 * t237 - t137 * t167 + (t154 * t237 + (-0.2e1 * t169 * t176 ^ 2 - t167) * t152) * qJD(4)) * t171) * t140;
	t1 = [t253 * t140 + t200 * t252, t117, 0, t115, 0, 0; t152 * t203 + (-t136 * t130 + (t118 * t152 + t121 * t134) * t131) * t126 + (t121 * t202 + (0.2e1 * t121 * t250 + (t134 * t140 - t134 - (-t122 * t140 * t216 + t224) * t155) * t131 * t138 + (-(-0.2e1 * t152 * t200 - t122) * t248 + (-(t122 + t218) * t155 + t253 * t152) * t131 * t140) * t139) * t126) * t155, t120 * t155 * t202 + (-(-t117 * t242 + (t122 * t245 - t136 * t139) * t129) * t248 + (t118 * t222 + t247) * t120 + (-t130 * t234 - (-t129 * t243 + t138 * t213 + t240) * t248) * t227) * t126 + (t203 * t234 + ((-t130 * t229 - (t206 * qJD(2) - t122) * t220) * t177 + (t130 * t205 + (t181 * t118 - (-t117 + t231) * t244 - (t206 * t122 - qJD(2)) * t241) * t131) * t175) * t126) * t174, 0, (t119 * t248 - t130 * t156) * t225 + (t119 * t247 + t135 * t130 + (t119 * t222 - t156 * t131) * t118 - (t176 * t230 - t175 * t228 - t115 * t152 - t125 * t136 + (-t125 * t235 - t154) * t122) * t131 * t241 - (-t137 + (-t115 * t174 - t122 * t176) * t175 - t189 * t125) * t220) * t126, 0, 0; t193 * t175 * t223 + (-t193 * t230 + ((qJD(1) * t149 + 0.2e1 * t154 * t246) * t181 + (-t251 * t135 - t137 * t181 + t154 * t205) * t150) * t175) * t144, (t149 * t177 * t181 + t176 * t217) * t223 + (0.2e1 * t176 * t199 + t191 * t149 + ((t135 * t181 - 0.2e1 * t176 * t210) * t177 + (t180 * t228 + 0.2e1 * t181 * t198) * t170) * t150) * t144, 0, t150 * t201 * t254 + (t201 * t246 + (t134 * t234 + (t175 * t205 - t177 * t229) * t155) * t150) * t144, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end