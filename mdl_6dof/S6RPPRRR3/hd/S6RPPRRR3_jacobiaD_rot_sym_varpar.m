% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR3
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
%   Wie in S6RPPRRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:50
	% DurationCPUTime: 0.96s
	% Computational Cost: add. (1787->89), mult. (2519->202), div. (480->12), fcn. (2968->9), ass. (0->95)
	t121 = qJ(1) + pkin(10);
	t120 = cos(t121);
	t193 = 0.2e1 * t120;
	t119 = sin(t121);
	t127 = sin(qJ(5));
	t128 = sin(qJ(4));
	t129 = cos(qJ(5));
	t171 = t128 * t129;
	t107 = t119 * t171 + t120 * t127;
	t104 = 0.1e1 / t107 ^ 2;
	t172 = t127 * t128;
	t106 = t119 * t172 - t120 * t129;
	t178 = t106 * t129;
	t103 = 0.1e1 / t107;
	t180 = t103 * t127;
	t141 = t104 * t178 - t180;
	t102 = t106 ^ 2;
	t101 = t102 * t104 + 0.1e1;
	t97 = 0.1e1 / t101;
	t192 = t141 * t97;
	t167 = qJD(4) * t120;
	t118 = t120 ^ 2;
	t123 = 0.1e1 / t128 ^ 2;
	t130 = cos(qJ(4));
	t126 = t130 ^ 2;
	t173 = t123 * t126;
	t115 = t118 * t173 + 0.1e1;
	t113 = 0.1e1 / t115;
	t122 = 0.1e1 / t128;
	t156 = t123 * t167;
	t168 = qJD(1) * t130;
	t158 = t119 * t168;
	t166 = qJD(4) * t128;
	t87 = ((t120 * t166 + t158) * t122 + t126 * t156) * t113;
	t150 = -t87 + t167;
	t151 = -t120 * t87 + qJD(4);
	t175 = t120 * t130;
	t112 = atan2(-t175, t128);
	t110 = sin(t112);
	t111 = cos(t112);
	t99 = -t110 * t175 + t111 * t128;
	t94 = 0.1e1 / t99;
	t95 = 0.1e1 / t99 ^ 2;
	t191 = t113 - 0.1e1;
	t117 = t119 ^ 2;
	t169 = qJD(1) * t120;
	t145 = t119 * t126 * t169;
	t160 = t95 * t166;
	t176 = t111 * t130;
	t82 = t151 * t176 + (t150 * t128 + t158) * t110;
	t189 = t82 * t94 * t95;
	t92 = t117 * t126 * t95 + 0.1e1;
	t190 = (t95 * t145 + (-t126 * t189 - t130 * t160) * t117) / t92 ^ 2;
	t179 = t104 * t106;
	t149 = qJD(5) * t128 + qJD(1);
	t165 = qJD(4) * t130;
	t138 = -t149 * t127 + t129 * t165;
	t148 = qJD(1) * t128 + qJD(5);
	t142 = t120 * t148;
	t89 = t138 * t119 + t129 * t142;
	t185 = t103 * t104 * t89;
	t139 = t127 * t165 + t149 * t129;
	t88 = t139 * t119 + t127 * t142;
	t188 = (-t102 * t185 + t88 * t179) / t101 ^ 2;
	t90 = 0.1e1 / t92;
	t187 = t90 * t95;
	t186 = t94 * t90;
	t125 = t130 * t126;
	t140 = qJD(4) * (-t123 * t125 - t130) * t122;
	t184 = 0.1e1 / t115 ^ 2 * (t118 * t140 - t123 * t145);
	t183 = t119 * t95;
	t177 = t110 * t128;
	t174 = t122 * t126;
	t170 = qJD(1) * t119;
	t164 = -0.2e1 * t189;
	t163 = 0.2e1 * t188;
	t162 = t94 * t190;
	t161 = t130 * t184;
	t159 = t113 * t174;
	t157 = t120 * t168;
	t155 = -0.2e1 * t95 * t190;
	t154 = 0.1e1 + t173;
	t153 = 0.2e1 * t106 * t185;
	t152 = t184 * t193;
	t147 = t120 * t159;
	t146 = t191 * t130 * t110;
	t144 = t154 * t119;
	t143 = t119 * t148;
	t109 = -t119 * t127 + t120 * t171;
	t108 = t119 * t129 + t120 * t172;
	t100 = t154 * t120 * t113;
	t86 = (-t111 * t147 - t146) * t119;
	t85 = t120 * t177 + t176 + (-t111 * t175 - t177) * t100;
	t83 = -t154 * t152 + (-qJD(1) * t144 + t140 * t193) * t113;
	t1 = [-0.2e1 * t119 * t122 * t161 + (-qJD(4) * t144 + t122 * t157) * t113, 0, 0, t83, 0, 0; (t166 * t186 + (0.2e1 * t162 + (qJD(1) * t86 + t82) * t187) * t130) * t120 + (t86 * t155 * t130 + (-t86 * t160 + (t86 * t164 + ((t87 * t147 + t191 * t166 + 0.2e1 * t161) * t110 + (t152 * t174 + t87 * t130 + (t125 * t156 + (-t87 + 0.2e1 * t167) * t130) * t113) * t111) * t183) * t130 + (t94 + ((t117 - t118) * t111 * t159 - t120 * t146) * t95) * t168) * t90) * t119, 0, 0, (t169 * t186 + (-0.2e1 * t162 + (-qJD(4) * t85 - t82) * t187) * t119) * t128 + (t85 * t119 * t155 + (t119 * qJD(4) * t94 + (t119 * t164 + t95 * t169) * t85 + (((t100 * t170 - t120 * t83) * t111 + (-t151 * t100 + t150) * t110) * t130 + ((-t83 - t170) * t110 + (t150 * t100 - t151) * t111) * t128) * t183) * t90) * t130, 0, 0; (-t103 * t108 + t109 * t179) * t163 + (t109 * t153 - t143 * t180 + t139 * t103 * t120 + (-t138 * t106 * t120 - t108 * t89 - t109 * t88 + t143 * t178) * t104) * t97, 0, 0, -t157 * t192 + (t166 * t192 + (t141 * t163 + ((qJD(5) * t103 + t153) * t129 + (-t129 * t88 + (qJD(5) * t106 - t89) * t127) * t104) * t97) * t130) * t119, -0.2e1 * t188 + 0.2e1 * (t104 * t88 * t97 + (-t104 * t188 - t97 * t185) * t106) * t106, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:50
	% DurationCPUTime: 1.05s
	% Computational Cost: add. (2446->95), mult. (2734->206), div. (498->12), fcn. (3199->9), ass. (0->97)
	t155 = sin(qJ(4));
	t150 = 0.1e1 / t155 ^ 2;
	t156 = cos(qJ(4));
	t153 = t156 ^ 2;
	t193 = t150 * t153;
	t148 = qJ(1) + pkin(10);
	t144 = cos(t148);
	t217 = 0.2e1 * t144;
	t216 = t156 * t193;
	t197 = t144 * t156;
	t136 = atan2(-t197, t155);
	t134 = sin(t136);
	t135 = cos(t136);
	t124 = -t134 * t197 + t135 * t155;
	t121 = 0.1e1 / t124;
	t143 = sin(t148);
	t154 = qJ(5) + qJ(6);
	t145 = sin(t154);
	t146 = cos(t154);
	t195 = t146 * t155;
	t131 = t143 * t195 + t144 * t145;
	t127 = 0.1e1 / t131;
	t149 = 0.1e1 / t155;
	t122 = 0.1e1 / t124 ^ 2;
	t128 = 0.1e1 / t131 ^ 2;
	t142 = t144 ^ 2;
	t139 = t142 * t193 + 0.1e1;
	t137 = 0.1e1 / t139;
	t215 = t137 - 0.1e1;
	t141 = t143 ^ 2;
	t199 = t141 * t153;
	t116 = t122 * t199 + 0.1e1;
	t191 = qJD(1) * t144;
	t171 = t143 * t153 * t191;
	t187 = qJD(4) * t156;
	t190 = qJD(1) * t156;
	t180 = t143 * t190;
	t188 = qJD(4) * t155;
	t189 = qJD(4) * t144;
	t113 = ((t144 * t188 + t180) * t149 + t189 * t193) * t137;
	t200 = t135 * t156;
	t107 = (-t113 * t144 + qJD(4)) * t200 + (t180 + (-t113 + t189) * t155) * t134;
	t212 = t107 * t121 * t122;
	t214 = (-t199 * t212 + (-t141 * t155 * t187 + t171) * t122) / t116 ^ 2;
	t147 = qJD(5) + qJD(6);
	t175 = t147 * t155 + qJD(1);
	t165 = t145 * t187 + t175 * t146;
	t174 = qJD(1) * t155 + t147;
	t168 = t144 * t174;
	t111 = t165 * t143 + t145 * t168;
	t196 = t145 * t155;
	t130 = t143 * t196 - t144 * t146;
	t126 = t130 ^ 2;
	t119 = t126 * t128 + 0.1e1;
	t203 = t128 * t130;
	t164 = -t175 * t145 + t146 * t187;
	t112 = t164 * t143 + t146 * t168;
	t211 = t112 * t127 * t128;
	t213 = (t111 * t203 - t126 * t211) / t119 ^ 2;
	t210 = t113 * t134;
	t209 = t113 * t156;
	t166 = qJD(4) * (-t156 - t216) * t149;
	t208 = (t142 * t166 - t150 * t171) / t139 ^ 2;
	t207 = t122 * t143;
	t206 = t122 * t156;
	t178 = 0.1e1 + t193;
	t125 = t178 * t144 * t137;
	t205 = t125 * t144;
	t204 = t127 * t145;
	t202 = t130 * t146;
	t201 = t134 * t155;
	t198 = t143 * t156;
	t194 = t149 * t153;
	t192 = qJD(1) * t143;
	t186 = 0.2e1 * t213;
	t185 = -0.2e1 * t212;
	t184 = t156 * t214;
	t183 = t156 * t208;
	t182 = t122 * t198;
	t181 = t137 * t194;
	t179 = t144 * t190;
	t177 = 0.2e1 * t130 * t211;
	t176 = t208 * t217;
	t173 = t144 * t181;
	t172 = t215 * t156 * t134;
	t170 = t178 * t143;
	t169 = t143 * t174;
	t167 = t128 * t202 - t204;
	t133 = -t143 * t145 + t144 * t195;
	t132 = t143 * t146 + t144 * t196;
	t117 = 0.1e1 / t119;
	t114 = 0.1e1 / t116;
	t110 = (-t135 * t173 - t172) * t143;
	t109 = t144 * t201 + t200 + (-t135 * t197 - t201) * t125;
	t108 = -t178 * t176 + (-qJD(1) * t170 + t166 * t217) * t137;
	t104 = -0.2e1 * t213 + 0.2e1 * (t111 * t117 * t128 + (-t117 * t211 - t128 * t213) * t130) * t130;
	t1 = [-0.2e1 * t143 * t149 * t183 + (-qJD(4) * t170 + t149 * t179) * t137, 0, 0, t108, 0, 0; (0.2e1 * t121 * t184 + (t121 * t188 + (qJD(1) * t110 + t107) * t206) * t114) * t144 + (-0.2e1 * t122 * t184 * t110 + (((t113 * t173 + t215 * t188 + 0.2e1 * t183) * t134 + (t176 * t194 + t209 + (-t209 + (0.2e1 * t156 + t216) * t189) * t137) * t135) * t182 + (-t122 * t188 + t156 * t185) * t110 + (t121 + ((t141 - t142) * t135 * t181 - t144 * t172) * t122) * t190) * t114) * t143, 0, 0, 0.2e1 * (-t109 * t206 - t121 * t155) * t143 * t214 + ((t121 * t191 + (-qJD(4) * t109 - t107) * t207) * t155 + (t143 * qJD(4) * t121 + (-t108 * t135 * t144 + t134 * t189 + t205 * t210 - t210 + (-qJD(4) * t134 + t135 * t192) * t125) * t182 + (t122 * t191 + t143 * t185) * t109 + ((-t108 - t192) * t134 + ((-0.1e1 + t205) * qJD(4) + (-t125 + t144) * t113) * t135) * t155 * t207) * t156) * t114, 0, 0; (-t127 * t132 + t133 * t203) * t186 + (t133 * t177 - t169 * t204 + t165 * t127 * t144 + (-t164 * t130 * t144 - t133 * t111 - t132 * t112 + t169 * t202) * t128) * t117, 0, 0, t167 * t186 * t198 + (-t167 * t179 + (t167 * t188 + ((t127 * t147 + t177) * t146 + (-t111 * t146 + (t130 * t147 - t112) * t145) * t128) * t156) * t143) * t117, t104, t104;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end