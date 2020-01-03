% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRP10
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
%   Wie in S5RPRRP10_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRRP10_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP10_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:52:35
	% EndTime: 2019-12-31 18:52:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:52:35
	% EndTime: 2019-12-31 18:52:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:52:35
	% EndTime: 2019-12-31 18:52:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:52:35
	% EndTime: 2019-12-31 18:52:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:52:36
	% EndTime: 2019-12-31 18:52:36
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (2270->94), mult. (2519->210), div. (480->12), fcn. (2968->9), ass. (0->93)
	t136 = sin(qJ(1));
	t133 = t136 ^ 2;
	t132 = pkin(8) + qJ(3);
	t130 = sin(t132);
	t126 = t130 ^ 2;
	t131 = cos(t132);
	t128 = 0.1e1 / t131 ^ 2;
	t185 = t126 * t128;
	t121 = t133 * t185 + 0.1e1;
	t125 = t130 * t126;
	t127 = 0.1e1 / t131;
	t182 = t127 * t130;
	t146 = qJD(3) * (t125 * t127 * t128 + t182);
	t138 = cos(qJ(1));
	t174 = qJD(1) * t138;
	t183 = t126 * t136;
	t151 = t174 * t183;
	t191 = (t128 * t151 + t133 * t146) / t121 ^ 2;
	t201 = -0.2e1 * t191;
	t158 = 0.1e1 + t185;
	t200 = t136 * t158;
	t137 = cos(qJ(4));
	t176 = t137 * t138;
	t135 = sin(qJ(4));
	t178 = t136 * t135;
	t117 = t131 * t176 + t178;
	t179 = t136 * t130;
	t118 = atan2(-t179, -t131);
	t113 = cos(t118);
	t112 = sin(t118);
	t164 = t112 * t179;
	t102 = -t113 * t131 - t164;
	t99 = 0.1e1 / t102;
	t109 = 0.1e1 / t117;
	t100 = 0.1e1 / t102 ^ 2;
	t110 = 0.1e1 / t117 ^ 2;
	t119 = 0.1e1 / t121;
	t199 = t119 - 0.1e1;
	t134 = t138 ^ 2;
	t173 = qJD(3) * t131;
	t184 = t126 * t134;
	t172 = qJD(3) * t136;
	t160 = t128 * t172;
	t161 = t130 * t174;
	t93 = (-(-t131 * t172 - t161) * t127 + t126 * t160) * t119;
	t155 = t93 - t172;
	t156 = -t136 * t93 + qJD(3);
	t187 = t113 * t130;
	t88 = t156 * t187 + (t155 * t131 - t161) * t112;
	t196 = t99 * t100 * t88;
	t96 = t100 * t184 + 0.1e1;
	t198 = (-t184 * t196 + (t130 * t134 * t173 - t151) * t100) / t96 ^ 2;
	t94 = 0.1e1 / t96;
	t197 = t100 * t94;
	t177 = t136 * t137;
	t180 = t135 * t138;
	t116 = t131 * t180 - t177;
	t108 = t116 ^ 2;
	t107 = t108 * t110 + 0.1e1;
	t189 = t110 * t116;
	t153 = -qJD(1) * t131 + qJD(4);
	t154 = qJD(4) * t131 - qJD(1);
	t171 = qJD(3) * t138;
	t159 = t130 * t171;
	t98 = -t154 * t180 + (t153 * t136 - t159) * t137;
	t194 = t109 * t110 * t98;
	t147 = t131 * t178 + t176;
	t97 = t147 * qJD(1) - t117 * qJD(4) + t135 * t159;
	t195 = 0.1e1 / t107 ^ 2 * (-t108 * t194 - t97 * t189);
	t190 = t109 * t135;
	t188 = t112 * t131;
	t186 = t116 * t137;
	t181 = t130 * t138;
	t175 = qJD(1) * t136;
	t170 = 0.2e1 * t198;
	t169 = 0.2e1 * t196;
	t168 = -0.2e1 * t195;
	t167 = t99 * t198;
	t166 = t116 * t194;
	t165 = t94 * t173;
	t163 = t119 * t126 * t127;
	t157 = t127 * t201;
	t152 = t136 * t163;
	t150 = t158 * t138;
	t149 = t153 * t138;
	t148 = t110 * t186 - t190;
	t115 = -t131 * t177 + t180;
	t105 = 0.1e1 / t107;
	t104 = t119 * t200;
	t92 = (t199 * t130 * t112 - t113 * t152) * t138;
	t90 = -t136 * t188 + t187 + (-t113 * t179 + t188) * t104;
	t89 = t200 * t201 + (qJD(1) * t150 + 0.2e1 * t136 * t146) * t119;
	t1 = [t157 * t181 + (qJD(3) * t150 - t175 * t182) * t119, 0, t89, 0, 0; (-t99 * t165 + (0.2e1 * t167 + (qJD(1) * t92 + t88) * t197) * t130) * t136 + ((-t92 * t165 + (t92 * t170 + ((0.2e1 * t130 * t191 - t93 * t152 - t199 * t173) * t112 + (t157 * t183 + t130 * t93 + (t125 * t160 - (t93 - 0.2e1 * t172) * t130) * t119) * t113) * t94 * t138) * t130) * t100 + (t92 * t169 + (-t99 + ((-t133 + t134) * t113 * t163 + t199 * t164) * t100) * qJD(1)) * t130 * t94) * t138, 0, (-t99 * t94 * t175 + (-0.2e1 * t167 + (-qJD(3) * t90 - t88) * t197) * t138) * t131 + (((-qJD(3) * t99 + t90 * t169) * t138 + (t90 * t175 + (-(-t104 * t174 - t136 * t89) * t113 - ((t104 * t136 - 0.1e1) * t93 + (-t104 + t136) * qJD(3)) * t112) * t181) * t100) * t94 + (t90 * t170 - ((t89 - t174) * t112 + (t155 * t104 + t156) * t113) * t94 * t131) * t100 * t138) * t130, 0, 0; 0.2e1 * (t109 * t147 + t115 * t189) * t195 + (0.2e1 * t115 * t166 - t154 * t109 * t177 + (t130 * t172 + t149) * t190 + (t147 * t98 + t115 * t97 - t149 * t186 - (qJD(3) * t130 * t137 + t154 * t135) * t116 * t136) * t110) * t105, 0, t148 * t168 * t181 + (t148 * t131 * t171 + (-t148 * t175 + ((-qJD(4) * t109 - 0.2e1 * t166) * t137 + (-t137 * t97 + (-qJD(4) * t116 + t98) * t135) * t110) * t138) * t130) * t105, t168 + 0.2e1 * (-t105 * t110 * t97 + (-t105 * t194 - t110 * t195) * t116) * t116, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:52:36
	% EndTime: 2019-12-31 18:52:36
	% DurationCPUTime: 0.73s
	% Computational Cost: add. (2270->94), mult. (2519->210), div. (480->12), fcn. (2968->9), ass. (0->93)
	t140 = sin(qJ(1));
	t137 = t140 ^ 2;
	t136 = pkin(8) + qJ(3);
	t134 = sin(t136);
	t130 = t134 ^ 2;
	t135 = cos(t136);
	t132 = 0.1e1 / t135 ^ 2;
	t189 = t130 * t132;
	t125 = t137 * t189 + 0.1e1;
	t129 = t134 * t130;
	t131 = 0.1e1 / t135;
	t186 = t131 * t134;
	t150 = qJD(3) * (t129 * t131 * t132 + t186);
	t142 = cos(qJ(1));
	t178 = qJD(1) * t142;
	t187 = t130 * t140;
	t155 = t178 * t187;
	t195 = (t132 * t155 + t137 * t150) / t125 ^ 2;
	t205 = -0.2e1 * t195;
	t162 = 0.1e1 + t189;
	t204 = t140 * t162;
	t141 = cos(qJ(4));
	t180 = t141 * t142;
	t139 = sin(qJ(4));
	t182 = t140 * t139;
	t121 = t135 * t180 + t182;
	t183 = t140 * t134;
	t122 = atan2(-t183, -t135);
	t117 = cos(t122);
	t116 = sin(t122);
	t168 = t116 * t183;
	t106 = -t117 * t135 - t168;
	t103 = 0.1e1 / t106;
	t113 = 0.1e1 / t121;
	t104 = 0.1e1 / t106 ^ 2;
	t114 = 0.1e1 / t121 ^ 2;
	t123 = 0.1e1 / t125;
	t203 = t123 - 0.1e1;
	t138 = t142 ^ 2;
	t188 = t130 * t138;
	t100 = t104 * t188 + 0.1e1;
	t177 = qJD(3) * t135;
	t176 = qJD(3) * t140;
	t164 = t132 * t176;
	t165 = t134 * t178;
	t97 = (-(-t135 * t176 - t165) * t131 + t130 * t164) * t123;
	t159 = t97 - t176;
	t160 = -t140 * t97 + qJD(3);
	t191 = t117 * t134;
	t92 = t160 * t191 + (t159 * t135 - t165) * t116;
	t200 = t103 * t104 * t92;
	t202 = (-t188 * t200 + (t134 * t138 * t177 - t155) * t104) / t100 ^ 2;
	t98 = 0.1e1 / t100;
	t201 = t104 * t98;
	t151 = t135 * t182 + t180;
	t175 = qJD(3) * t142;
	t163 = t134 * t175;
	t101 = t151 * qJD(1) - qJD(4) * t121 + t139 * t163;
	t181 = t140 * t141;
	t184 = t139 * t142;
	t120 = t135 * t184 - t181;
	t112 = t120 ^ 2;
	t111 = t112 * t114 + 0.1e1;
	t193 = t114 * t120;
	t157 = -qJD(1) * t135 + qJD(4);
	t158 = qJD(4) * t135 - qJD(1);
	t102 = -t158 * t184 + (t157 * t140 - t163) * t141;
	t197 = t102 * t113 * t114;
	t199 = 0.1e1 / t111 ^ 2 * (-t101 * t193 - t112 * t197);
	t194 = t113 * t139;
	t192 = t116 * t135;
	t190 = t120 * t141;
	t185 = t134 * t142;
	t179 = qJD(1) * t140;
	t174 = 0.2e1 * t202;
	t173 = 0.2e1 * t200;
	t172 = -0.2e1 * t199;
	t171 = t103 * t202;
	t170 = t98 * t177;
	t169 = t120 * t197;
	t167 = t123 * t130 * t131;
	t161 = t131 * t205;
	t156 = t140 * t167;
	t154 = t162 * t142;
	t153 = t157 * t142;
	t152 = t114 * t190 - t194;
	t119 = -t135 * t181 + t184;
	t109 = 0.1e1 / t111;
	t108 = t123 * t204;
	t96 = (t203 * t134 * t116 - t117 * t156) * t142;
	t94 = -t140 * t192 + t191 + (-t117 * t183 + t192) * t108;
	t93 = t204 * t205 + (qJD(1) * t154 + 0.2e1 * t140 * t150) * t123;
	t1 = [t161 * t185 + (qJD(3) * t154 - t179 * t186) * t123, 0, t93, 0, 0; (-t103 * t170 + (0.2e1 * t171 + (qJD(1) * t96 + t92) * t201) * t134) * t140 + ((-t96 * t170 + (t96 * t174 + ((0.2e1 * t134 * t195 - t97 * t156 - t203 * t177) * t116 + (t161 * t187 + t134 * t97 + (t129 * t164 - (t97 - 0.2e1 * t176) * t134) * t123) * t117) * t98 * t142) * t134) * t104 + (t96 * t173 + (-t103 + ((-t137 + t138) * t117 * t167 + t203 * t168) * t104) * qJD(1)) * t134 * t98) * t142, 0, (-t103 * t98 * t179 + (-0.2e1 * t171 + (-qJD(3) * t94 - t92) * t201) * t142) * t135 + (((-qJD(3) * t103 + t94 * t173) * t142 + (t94 * t179 + (-(-t108 * t178 - t140 * t93) * t117 - ((t108 * t140 - 0.1e1) * t97 + (-t108 + t140) * qJD(3)) * t116) * t185) * t104) * t98 + (t94 * t174 - ((t93 - t178) * t116 + (t159 * t108 + t160) * t117) * t98 * t135) * t104 * t142) * t134, 0, 0; 0.2e1 * (t113 * t151 + t119 * t193) * t199 + (0.2e1 * t119 * t169 - t158 * t113 * t181 + (t134 * t176 + t153) * t194 + (t119 * t101 + t151 * t102 - t153 * t190 - (qJD(3) * t134 * t141 + t158 * t139) * t120 * t140) * t114) * t109, 0, t152 * t172 * t185 + (t152 * t135 * t175 + (-t152 * t179 + ((-qJD(4) * t113 - 0.2e1 * t169) * t141 + (-t101 * t141 + (-qJD(4) * t120 + t102) * t139) * t114) * t142) * t134) * t109, t172 + 0.2e1 * (-t101 * t109 * t114 + (-t109 * t197 - t114 * t199) * t120) * t120, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end