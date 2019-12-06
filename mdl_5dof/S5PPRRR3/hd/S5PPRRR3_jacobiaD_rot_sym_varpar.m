% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPRRR3
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
%   Wie in S5PPRRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PPRRR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:17:26
	% EndTime: 2019-12-05 15:17:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:17:26
	% EndTime: 2019-12-05 15:17:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:17:26
	% EndTime: 2019-12-05 15:17:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:17:26
	% EndTime: 2019-12-05 15:17:27
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t43 = sin(pkin(8));
	t46 = sin(qJ(3));
	t47 = cos(qJ(3));
	t50 = cos(pkin(8)) * cos(pkin(9));
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
	% StartTime: 2019-12-05 15:17:27
	% EndTime: 2019-12-05 15:17:27
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (682->54), mult. (2271->129), div. (423->14), fcn. (2956->11), ass. (0->67)
	t143 = cos(pkin(9));
	t144 = cos(pkin(8));
	t148 = cos(qJ(3));
	t162 = t144 * t148;
	t142 = sin(pkin(8));
	t146 = sin(qJ(3));
	t165 = t142 * t146;
	t131 = t143 * t162 + t165;
	t145 = sin(qJ(4));
	t147 = cos(qJ(4));
	t141 = sin(pkin(9));
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
	% StartTime: 2019-12-05 15:17:27
	% EndTime: 2019-12-05 15:17:27
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (1083->57), mult. (2577->133), div. (441->14), fcn. (3313->11), ass. (0->73)
	t179 = cos(pkin(9));
	t180 = cos(pkin(8));
	t182 = cos(qJ(3));
	t198 = t180 * t182;
	t178 = sin(pkin(8));
	t181 = sin(qJ(3));
	t201 = t178 * t181;
	t159 = t179 * t201 + t198;
	t177 = sin(pkin(9));
	t202 = t177 * t181;
	t151 = atan2(-t159, t202);
	t147 = sin(t151);
	t148 = cos(t151);
	t134 = -t147 * t159 + t148 * t202;
	t131 = 0.1e1 / t134;
	t163 = t179 * t198 + t201;
	t176 = qJ(4) + qJ(5);
	t168 = sin(t176);
	t169 = cos(t176);
	t203 = t177 * t180;
	t146 = t163 * t169 + t168 * t203;
	t142 = 0.1e1 / t146;
	t173 = 0.1e1 / t181;
	t132 = 0.1e1 / t134 ^ 2;
	t143 = 0.1e1 / t146 ^ 2;
	t174 = 0.1e1 / t181 ^ 2;
	t200 = t178 * t182;
	t193 = t179 * t200;
	t197 = qJD(3) * t181;
	t154 = qJD(3) * t193 - t180 * t197;
	t204 = t174 * t182;
	t194 = t159 * t204;
	t157 = t159 ^ 2;
	t171 = 0.1e1 / t177 ^ 2;
	t152 = t157 * t171 * t174 + 0.1e1;
	t149 = 0.1e1 / t152;
	t170 = 0.1e1 / t177;
	t206 = t149 * t170;
	t126 = (qJD(3) * t194 - t154 * t173) * t206;
	t190 = -t147 * t202 - t148 * t159;
	t195 = t148 * t177 * t182;
	t123 = qJD(3) * t195 + t190 * t126 - t147 * t154;
	t211 = t123 * t131 * t132;
	t199 = t180 * t181;
	t162 = t179 * t199 - t200;
	t210 = t132 * t162;
	t145 = t163 * t168 - t169 * t203;
	t141 = t145 ^ 2;
	t138 = t141 * t143 + 0.1e1;
	t155 = t162 * qJD(3);
	t172 = qJD(4) + qJD(5);
	t192 = t172 * t203 - t155;
	t205 = t163 * t172;
	t139 = t192 * t168 + t169 * t205;
	t207 = t143 * t145;
	t140 = -t168 * t205 + t192 * t169;
	t208 = t140 * t142 * t143;
	t209 = 0.1e1 / t138 ^ 2 * (t139 * t207 - t141 * t208);
	t196 = -0.2e1 * t209;
	t191 = -t142 * t168 + t169 * t207;
	t161 = t193 - t199;
	t189 = -t161 * t173 + t194;
	t175 = t173 * t174;
	t158 = t162 ^ 2;
	t156 = t163 * qJD(3);
	t153 = t159 * qJD(3);
	t135 = 0.1e1 / t138;
	t130 = t132 * t158 + 0.1e1;
	t127 = t189 * t206;
	t124 = t190 * t127 - t147 * t161 + t195;
	t122 = (-0.2e1 * t189 / t152 ^ 2 * (-qJD(3) * t157 * t175 * t182 + t154 * t159 * t174) * t171 + (t154 * t204 + t153 * t173 + (t161 * t204 + (-0.2e1 * t175 * t182 ^ 2 - t173) * t159) * qJD(3)) * t149) * t170;
	t120 = t196 + 0.2e1 * (t135 * t139 * t143 + (-t135 * t208 - t143 * t209) * t145) * t145;
	t1 = [0, 0, t122, 0, 0; 0, 0, 0.2e1 * (t124 * t210 - t131 * t163) / t130 ^ 2 * (t156 * t210 - t158 * t211) + (-t155 * t131 + (-t163 * t123 - t124 * t156) * t132 + (0.2e1 * t124 * t211 + (-(-t177 * t197 - t122 * t159 - t127 * t154 + (-t127 * t202 - t161) * t126) * t148 - (t126 * t127 * t159 + t153 + (-t122 * t181 + (-qJD(3) * t127 - t126) * t182) * t177) * t147) * t132) * t162) / t130, 0, 0; 0, 0, t191 * t162 * t196 + (t191 * t156 + ((-t142 * t172 - 0.2e1 * t145 * t208) * t169 + (t139 * t169 + (-t145 * t172 + t140) * t168) * t143) * t162) * t135, t120, t120;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end