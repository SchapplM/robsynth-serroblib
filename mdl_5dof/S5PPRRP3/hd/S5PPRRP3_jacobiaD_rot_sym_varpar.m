% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPRRP3
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
%   Wie in S5PPRRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PPRRP3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:38
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t43 = sin(pkin(7));
	t46 = sin(qJ(3));
	t47 = cos(qJ(3));
	t50 = cos(pkin(7)) * cos(pkin(8));
	t41 = t43 * t46 + t47 * t50;
	t38 = 0.1e1 / t41 ^ 2;
	t54 = qJD(3) * t38;
	t40 = -t43 * t47 + t46 * t50;
	t37 = t40 ^ 2;
	t34 = t37 * t38 + 0.1e1;
	t51 = t41 * t54;
	t52 = t40 / t41 * t54;
	t53 = (t37 * t52 + t40 * t51) / t34 ^ 2;
	t32 = 0.1e1 / t34;
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, -0.2e1 * t53 + 0.2e1 * (t32 * t51 + (t32 * t52 - t38 * t53) * t40) * t40, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:39
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (682->54), mult. (2271->129), div. (423->14), fcn. (2956->11), ass. (0->67)
	t143 = cos(pkin(8));
	t144 = cos(pkin(7));
	t148 = cos(qJ(3));
	t162 = t144 * t148;
	t142 = sin(pkin(7));
	t146 = sin(qJ(3));
	t165 = t142 * t146;
	t131 = t143 * t162 + t165;
	t145 = sin(qJ(4));
	t147 = cos(qJ(4));
	t141 = sin(pkin(8));
	t167 = t141 * t144;
	t155 = -t131 * t145 + t147 * t167;
	t176 = t155 * qJD(4);
	t163 = t144 * t146;
	t164 = t142 * t148;
	t129 = t143 * t164 - t163;
	t127 = t143 * t165 + t162;
	t166 = t141 * t146;
	t119 = atan2(-t127, t166);
	t115 = sin(t119);
	t116 = cos(t119);
	t102 = -t115 * t127 + t116 * t166;
	t99 = 0.1e1 / t102;
	t114 = t131 * t147 + t145 * t167;
	t110 = 0.1e1 / t114;
	t138 = 0.1e1 / t146;
	t100 = 0.1e1 / t102 ^ 2;
	t111 = 0.1e1 / t114 ^ 2;
	t139 = 0.1e1 / t146 ^ 2;
	t122 = t129 * qJD(3);
	t161 = qJD(3) * t148;
	t168 = t139 * t148;
	t159 = t127 * t168;
	t125 = t127 ^ 2;
	t137 = 0.1e1 / t141 ^ 2;
	t120 = t125 * t137 * t139 + 0.1e1;
	t117 = 0.1e1 / t120;
	t136 = 0.1e1 / t141;
	t169 = t117 * t136;
	t94 = (qJD(3) * t159 - t122 * t138) * t169;
	t91 = (-t127 * t94 + t141 * t161) * t116 + (-t94 * t166 - t122) * t115;
	t175 = t99 * t100 * t91;
	t109 = t155 ^ 2;
	t106 = t109 * t111 + 0.1e1;
	t130 = t143 * t163 - t164;
	t123 = t130 * qJD(3);
	t107 = t114 * qJD(4) - t123 * t145;
	t170 = t111 * t155;
	t108 = -t123 * t147 + t176;
	t171 = t108 * t110 * t111;
	t174 = 0.1e1 / t106 ^ 2 * (-t107 * t170 - t109 * t171);
	t156 = -t129 * t138 + t159;
	t95 = t156 * t169;
	t173 = t127 * t95;
	t172 = t100 * t130;
	t160 = -0.2e1 * t174;
	t157 = -t110 * t145 - t147 * t170;
	t140 = t138 * t139;
	t126 = t130 ^ 2;
	t124 = t131 * qJD(3);
	t121 = t127 * qJD(3);
	t104 = 0.1e1 / t106;
	t98 = t100 * t126 + 0.1e1;
	t92 = (t141 * t148 - t173) * t116 + (-t95 * t166 - t129) * t115;
	t90 = (-0.2e1 * t156 / t120 ^ 2 * (t122 * t127 * t139 - t125 * t140 * t161) * t137 + (t122 * t168 + t121 * t138 + (t129 * t168 + (-0.2e1 * t140 * t148 ^ 2 - t138) * t127) * qJD(3)) * t117) * t136;
	t1 = [0, 0, t90, 0, 0; 0, 0, 0.2e1 * (-t131 * t99 + t92 * t172) / t98 ^ 2 * (t124 * t172 - t126 * t175) + (-t123 * t99 + (-t92 * t124 - t131 * t91) * t100 + (0.2e1 * t92 * t175 + (-(-t122 * t95 - t127 * t90 - t129 * t94 + (-t94 * t95 - qJD(3)) * t166) * t116 - (t94 * t173 + t121 + (-t146 * t90 + (-qJD(3) * t95 - t94) * t148) * t141) * t115) * t100) * t130) / t98, 0, 0; 0, 0, t157 * t130 * t160 + (t157 * t124 + ((-qJD(4) * t110 + 0.2e1 * t155 * t171) * t147 + (t107 * t147 + (t108 + t176) * t145) * t111) * t130) * t104, t160 - 0.2e1 * (t104 * t107 * t111 - (-t104 * t171 - t111 * t174) * t155) * t155, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:11:38
	% EndTime: 2019-12-05 15:11:39
	% DurationCPUTime: 0.78s
	% Computational Cost: add. (2420->86), mult. (7430->189), div. (524->12), fcn. (9601->11), ass. (0->90)
	t172 = cos(pkin(8));
	t177 = cos(qJ(3));
	t217 = sin(pkin(7));
	t189 = t217 * t177;
	t173 = cos(pkin(7));
	t175 = sin(qJ(3));
	t202 = t173 * t175;
	t163 = t172 * t189 - t202;
	t174 = sin(qJ(4));
	t176 = cos(qJ(4));
	t171 = sin(pkin(8));
	t191 = t171 * t217;
	t149 = t163 * t176 + t174 * t191;
	t190 = t217 * t175;
	t201 = t173 * t177;
	t162 = -t172 * t190 - t201;
	t154 = t162 * qJD(3);
	t129 = t149 * qJD(4) + t154 * t174;
	t184 = -t163 * t174 + t176 * t191;
	t142 = t184 ^ 2;
	t203 = t171 * t177;
	t166 = t172 * t176 + t174 * t203;
	t160 = 0.1e1 / t166 ^ 2;
	t140 = t142 * t160 + 0.1e1;
	t137 = 0.1e1 / t140;
	t167 = -t172 * t174 + t176 * t203;
	t204 = t171 * t175;
	t192 = qJD(3) * t204;
	t152 = t167 * qJD(4) - t174 * t192;
	t159 = 0.1e1 / t166;
	t210 = t184 * t160;
	t117 = (-t129 * t159 - t152 * t210) * t137;
	t141 = atan2(t184, t166);
	t133 = sin(t141);
	t134 = cos(t141);
	t188 = -t133 * t166 + t134 * t184;
	t114 = t188 * t117 - t133 * t129 + t134 * t152;
	t128 = t133 * t184 + t134 * t166;
	t125 = 0.1e1 / t128;
	t126 = 0.1e1 / t128 ^ 2;
	t222 = t114 * t125 * t126;
	t221 = 0.2e1 * t176;
	t165 = t172 * t201 + t190;
	t205 = t171 * t173;
	t186 = -t165 * t174 + t176 * t205;
	t220 = -0.2e1 * t186 * t222;
	t208 = t152 * t159 * t160;
	t219 = (-t129 * t210 - t142 * t208) / t140 ^ 2;
	t194 = t184 * t204;
	t185 = t159 * t162 - t160 * t194;
	t218 = t174 * t185;
	t164 = -t172 * t202 + t189;
	t151 = t165 * t176 + t174 * t205;
	t144 = 0.1e1 / t151;
	t145 = 0.1e1 / t151 ^ 2;
	t158 = t164 ^ 2;
	t211 = t145 * t158;
	t139 = 0.1e1 + t211;
	t156 = t164 * qJD(3);
	t132 = t186 * qJD(4) + t156 * t176;
	t214 = t132 * t144 * t145;
	t195 = t158 * t214;
	t157 = t165 * qJD(3);
	t207 = t157 * t164;
	t216 = (-t145 * t207 - t195) / t139 ^ 2;
	t215 = t126 * t186;
	t213 = t133 * t186;
	t212 = t134 * t186;
	t209 = t184 * t167;
	t206 = t164 * t174;
	t200 = qJD(3) * t177;
	t199 = qJD(4) * t176;
	t143 = t186 ^ 2;
	t123 = t126 * t143 + 0.1e1;
	t131 = t151 * qJD(4) + t156 * t174;
	t198 = 0.2e1 * (-t131 * t215 - t143 * t222) / t123 ^ 2;
	t196 = 0.2e1 * t186 * t164;
	t187 = -t149 * t159 - t160 * t209;
	t155 = t163 * qJD(3);
	t153 = -t166 * qJD(4) - t176 * t192;
	t135 = 0.1e1 / t139;
	t130 = t184 * qJD(4) + t154 * t176;
	t120 = 0.1e1 / t123;
	t119 = t137 * t218;
	t118 = t187 * t137;
	t116 = (-t133 * t162 - t134 * t204) * t174 - t188 * t119;
	t115 = t188 * t118 - t133 * t149 + t134 * t167;
	t113 = -0.2e1 * t187 * t219 + (0.2e1 * t208 * t209 - t130 * t159 + (t129 * t167 + t149 * t152 - t153 * t184) * t160) * t137;
	t111 = 0.2e1 * t218 * t219 + (-t185 * t199 + (-0.2e1 * t194 * t208 + t155 * t159 + (t152 * t162 + (-t129 * t175 + t184 * t200) * t171) * t160) * t174) * t137;
	t1 = [0, 0, t111, t113, 0; 0, 0, (-t116 * t215 - t125 * t206) * t198 + ((-t157 * t174 + t164 * t199) * t125 + t116 * t220 + (-t116 * t131 - t206 * t114 + (t111 * t184 + t119 * t129 + (-t174 * t200 - t175 * t199) * t171 + (t119 * t166 - t162 * t174) * t117) * t212 + (-t162 * t199 - t111 * t166 + t119 * t152 + t155 * t174 + (t119 * t184 + t174 * t204) * t117) * t213) * t126) * t120, (-t115 * t215 - t125 * t151) * t198 + (t115 * t220 + t132 * t125 + (-t151 * t114 - t115 * t131 + (t113 * t184 - t118 * t129 + t153 + (-t118 * t166 - t149) * t117) * t212 + (-t113 * t166 - t118 * t152 - t130 + (-t118 * t184 - t167) * t117) * t213) * t126) * t120, 0; 0, 0, 0.2e1 * (t144 * t165 + t176 * t211) * t216 + (t195 * t221 - t144 * t156 + (qJD(4) * t158 * t174 + t132 * t165 + t207 * t221) * t145) * t135, t145 * t196 * t216 + (t196 * t214 + (t131 * t164 + t157 * t186) * t145) * t135, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end