% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRP1
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
%   Wie in S6RPPRRP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:47
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRP1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:47:46
	% EndTime: 2019-10-09 23:47:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:47:46
	% EndTime: 2019-10-09 23:47:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:47:46
	% EndTime: 2019-10-09 23:47:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:47:46
	% EndTime: 2019-10-09 23:47:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:47:46
	% EndTime: 2019-10-09 23:47:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:47:47
	% EndTime: 2019-10-09 23:47:48
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (3244->94), mult. (2519->205), div. (480->12), fcn. (2968->9), ass. (0->93)
	t144 = qJ(1) + pkin(9);
	t140 = sin(t144);
	t134 = t140 ^ 2;
	t143 = pkin(10) + qJ(4);
	t139 = sin(t143);
	t133 = t139 ^ 2;
	t141 = cos(t143);
	t136 = 0.1e1 / t141 ^ 2;
	t193 = t133 * t136;
	t128 = t134 * t193 + 0.1e1;
	t132 = t139 * t133;
	t135 = 0.1e1 / t141;
	t190 = t135 * t139;
	t154 = qJD(4) * (t132 * t135 * t136 + t190);
	t142 = cos(t144);
	t181 = qJD(1) * t142;
	t191 = t133 * t140;
	t159 = t181 * t191;
	t197 = (t134 * t154 + t136 * t159) / t128 ^ 2;
	t207 = -0.2e1 * t197;
	t165 = 0.1e1 + t193;
	t206 = t140 * t165;
	t146 = cos(qJ(5));
	t184 = t142 * t146;
	t145 = sin(qJ(5));
	t187 = t140 * t145;
	t122 = t141 * t184 + t187;
	t162 = qJD(5) * t141 - qJD(1);
	t180 = qJD(4) * t139;
	t205 = t162 * t145 + t146 * t180;
	t188 = t140 * t139;
	t125 = atan2(-t188, -t141);
	t124 = cos(t125);
	t123 = sin(t125);
	t173 = t123 * t188;
	t109 = -t124 * t141 - t173;
	t106 = 0.1e1 / t109;
	t116 = 0.1e1 / t122;
	t107 = 0.1e1 / t109 ^ 2;
	t117 = 0.1e1 / t122 ^ 2;
	t126 = 0.1e1 / t128;
	t204 = t126 - 0.1e1;
	t172 = t126 * t133 * t135;
	t160 = t140 * t172;
	t99 = (t204 * t139 * t123 - t124 * t160) * t142;
	t203 = t107 * t99;
	t179 = qJD(4) * t140;
	t169 = t136 * t179;
	t170 = t139 * t181;
	t178 = qJD(4) * t141;
	t100 = (-(-t140 * t178 - t170) * t135 + t133 * t169) * t126;
	t194 = t124 * t139;
	t95 = (-t100 * t140 + qJD(4)) * t194 + (-t170 + (t100 - t179) * t141) * t123;
	t202 = t106 * t107 * t95;
	t155 = t141 * t187 + t184;
	t168 = t145 * t180;
	t104 = t155 * qJD(1) - t122 * qJD(5) + t142 * t168;
	t185 = t142 * t145;
	t186 = t140 * t146;
	t121 = t141 * t185 - t186;
	t115 = t121 ^ 2;
	t114 = t115 * t117 + 0.1e1;
	t196 = t117 * t121;
	t161 = -qJD(1) * t141 + qJD(5);
	t157 = t161 * t146;
	t105 = t140 * t157 - t205 * t142;
	t200 = t105 * t116 * t117;
	t201 = 0.1e1 / t114 ^ 2 * (-t104 * t196 - t115 * t200);
	t199 = t107 * t139;
	t198 = t107 * t142;
	t195 = t123 * t141;
	t138 = t142 ^ 2;
	t192 = t133 * t138;
	t189 = t139 * t142;
	t111 = t126 * t206;
	t183 = t111 - t140;
	t182 = qJD(1) * t140;
	t103 = t107 * t192 + 0.1e1;
	t177 = 0.2e1 / t103 ^ 2 * (-t192 * t202 + (t138 * t139 * t178 - t159) * t107);
	t176 = 0.2e1 * t202;
	t175 = -0.2e1 * t201;
	t174 = t121 * t200;
	t166 = t111 * t140 - 0.1e1;
	t164 = t139 * t177;
	t163 = t135 * t207;
	t158 = t165 * t142;
	t156 = -t116 * t145 + t146 * t196;
	t120 = -t141 * t186 + t185;
	t112 = 0.1e1 / t114;
	t101 = 0.1e1 / t103;
	t97 = -t140 * t195 + t194 + (-t124 * t188 + t195) * t111;
	t96 = t206 * t207 + (qJD(1) * t158 + 0.2e1 * t140 * t154) * t126;
	t1 = [t163 * t189 + (qJD(4) * t158 - t182 * t190) * t126, 0, 0, t96, 0, 0; (t106 * t164 + (-t106 * t178 + (qJD(1) * t99 + t95) * t199) * t101) * t140 + (t164 * t203 + (-t178 * t203 + (t99 * t176 + ((-t100 * t160 + 0.2e1 * t139 * t197 - t204 * t178) * t123 + (t163 * t191 + t100 * t139 + (t132 * t169 - (t100 - 0.2e1 * t179) * t139) * t126) * t124) * t198) * t139 + (-t106 + ((-t134 + t138) * t124 * t172 + t204 * t173) * t107) * t139 * qJD(1)) * t101) * t142, 0, 0, (-t106 * t141 + t97 * t199) * t142 * t177 + ((-t106 * t182 + (-qJD(4) * t97 - t95) * t198) * t141 + ((-qJD(4) * t106 + t97 * t176) * t142 + (t97 * t182 + (-(-t111 * t181 - t140 * t96) * t124 - (-t183 * qJD(4) + t166 * t100) * t123) * t189) * t107 - ((t96 - t181) * t123 + (-t166 * qJD(4) + t183 * t100) * t124) * t141 * t198) * t139) * t101, 0, 0; 0.2e1 * (t116 * t155 + t120 * t196) * t201 + (0.2e1 * t120 * t174 + (t120 * t104 + t155 * t105 + (-t205 * t140 - t142 * t157) * t121) * t117 + (t161 * t185 + (-t162 * t146 + t168) * t140) * t116) * t112, 0, 0, t156 * t175 * t189 + (t156 * t142 * t178 + (-t156 * t182 + ((-qJD(5) * t116 - 0.2e1 * t174) * t146 + (-t104 * t146 + (-qJD(5) * t121 + t105) * t145) * t117) * t142) * t139) * t112, t175 + 0.2e1 * (-t104 * t112 * t117 + (-t112 * t200 - t117 * t201) * t121) * t121, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:47:47
	% EndTime: 2019-10-09 23:47:48
	% DurationCPUTime: 1.11s
	% Computational Cost: add. (3244->96), mult. (2519->203), div. (480->12), fcn. (2968->9), ass. (0->95)
	t147 = pkin(10) + qJ(4);
	t143 = sin(t147);
	t137 = t143 ^ 2;
	t145 = cos(t147);
	t140 = 0.1e1 / t145 ^ 2;
	t196 = t137 * t140;
	t148 = qJ(1) + pkin(9);
	t144 = sin(t148);
	t214 = 0.2e1 * t144;
	t213 = t143 * t196;
	t146 = cos(t148);
	t150 = cos(qJ(5));
	t188 = t146 * t150;
	t149 = sin(qJ(5));
	t191 = t144 * t149;
	t126 = t145 * t188 + t191;
	t166 = qJD(5) * t145 - qJD(1);
	t185 = qJD(4) * t143;
	t212 = t166 * t149 + t150 * t185;
	t192 = t144 * t143;
	t129 = atan2(-t192, -t145);
	t128 = cos(t129);
	t127 = sin(t129);
	t176 = t127 * t192;
	t113 = -t128 * t145 - t176;
	t110 = 0.1e1 / t113;
	t120 = 0.1e1 / t126;
	t139 = 0.1e1 / t145;
	t111 = 0.1e1 / t113 ^ 2;
	t121 = 0.1e1 / t126 ^ 2;
	t211 = -0.2e1 * t143;
	t138 = t144 ^ 2;
	t132 = t138 * t196 + 0.1e1;
	t130 = 0.1e1 / t132;
	t210 = t130 - 0.1e1;
	t186 = qJD(1) * t146;
	t173 = t143 * t186;
	t183 = qJD(4) * t145;
	t184 = qJD(4) * t144;
	t104 = (-(-t144 * t183 - t173) * t139 + t184 * t196) * t130;
	t198 = t128 * t143;
	t99 = (-t104 * t144 + qJD(4)) * t198 + (-t173 + (t104 - t184) * t145) * t127;
	t209 = t110 * t111 * t99;
	t159 = t145 * t191 + t188;
	t172 = t149 * t185;
	t108 = t159 * qJD(1) - qJD(5) * t126 + t146 * t172;
	t189 = t146 * t149;
	t190 = t144 * t150;
	t125 = t145 * t189 - t190;
	t119 = t125 ^ 2;
	t118 = t119 * t121 + 0.1e1;
	t200 = t121 * t125;
	t165 = -qJD(1) * t145 + qJD(5);
	t161 = t165 * t150;
	t109 = t144 * t161 - t146 * t212;
	t205 = t109 * t120 * t121;
	t208 = (-t108 * t200 - t119 * t205) / t118 ^ 2;
	t207 = t104 * t127;
	t206 = t104 * t143;
	t204 = t111 * t143;
	t203 = t111 * t146;
	t194 = t139 * t143;
	t158 = qJD(4) * (t139 * t213 + t194);
	t163 = t137 * t144 * t186;
	t202 = (t138 * t158 + t140 * t163) / t132 ^ 2;
	t170 = 0.1e1 + t196;
	t115 = t170 * t144 * t130;
	t201 = t115 * t144;
	t199 = t127 * t145;
	t197 = t137 * t139;
	t142 = t146 ^ 2;
	t195 = t137 * t142;
	t193 = t143 * t146;
	t187 = qJD(1) * t144;
	t182 = qJD(4) * t146;
	t107 = t111 * t195 + 0.1e1;
	t181 = 0.2e1 / t107 ^ 2 * (-t195 * t209 + (t142 * t143 * t183 - t163) * t111);
	t180 = 0.2e1 * t209;
	t179 = -0.2e1 * t208;
	t178 = t125 * t205;
	t177 = t111 * t193;
	t175 = t130 * t197;
	t169 = t143 * t181;
	t168 = t202 * t211;
	t167 = t202 * t214;
	t164 = t144 * t175;
	t162 = t170 * t146;
	t160 = -t120 * t149 + t150 * t200;
	t124 = -t145 * t190 + t189;
	t116 = 0.1e1 / t118;
	t105 = 0.1e1 / t107;
	t103 = (t210 * t143 * t127 - t128 * t164) * t146;
	t101 = -t144 * t199 + t198 + (-t128 * t192 + t199) * t115;
	t100 = -t170 * t167 + (qJD(1) * t162 + t158 * t214) * t130;
	t1 = [t139 * t146 * t168 + (qJD(4) * t162 - t187 * t194) * t130, 0, 0, t100, 0, 0; (t110 * t169 + (-t110 * t183 + (qJD(1) * t103 + t99) * t204) * t105) * t144 + (t111 * t169 * t103 + (-((t104 * t164 + t210 * t183 + t168) * t127 + (t167 * t197 - t206 + (t206 + (t211 - t213) * t184) * t130) * t128) * t177 + (-t111 * t183 + t143 * t180) * t103 + (-t110 + ((-t138 + t142) * t128 * t175 + t210 * t176) * t111) * t143 * qJD(1)) * t105) * t146, 0, 0, (t101 * t204 - t110 * t145) * t146 * t181 + ((-t110 * t187 + (-qJD(4) * t101 - t99) * t203) * t145 + (-t110 * t182 - (-t100 * t128 * t144 + t127 * t184 + t201 * t207 - t207 + (-qJD(4) * t127 - t128 * t186) * t115) * t177 + (t111 * t187 + t146 * t180) * t101 - ((t100 - t186) * t127 + ((0.1e1 - t201) * qJD(4) + (t115 - t144) * t104) * t128) * t145 * t203) * t143) * t105, 0, 0; 0.2e1 * (t120 * t159 + t124 * t200) * t208 + (0.2e1 * t124 * t178 + (t124 * t108 + t159 * t109 + (-t144 * t212 - t146 * t161) * t125) * t121 + (t165 * t189 + (-t166 * t150 + t172) * t144) * t120) * t116, 0, 0, t160 * t179 * t193 + (t160 * t145 * t182 + (-t160 * t187 + ((-qJD(5) * t120 - 0.2e1 * t178) * t150 + (-t108 * t150 + (-qJD(5) * t125 + t109) * t149) * t121) * t146) * t143) * t116, t179 + 0.2e1 * (-t108 * t116 * t121 + (-t116 * t205 - t121 * t208) * t125) * t125, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end