% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S4RRRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4RRRP6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRRP6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_jacobiaD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:19:24
	% EndTime: 2019-12-31 17:19:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:19:24
	% EndTime: 2019-12-31 17:19:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:19:24
	% EndTime: 2019-12-31 17:19:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:19:24
	% EndTime: 2019-12-31 17:19:25
	% DurationCPUTime: 0.69s
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
	t1 = [t145 * t171 + (qJD(2) * t141 - t164 * t172) * t111, t81, 0, 0; (-t161 * t184 + (0.2e1 * t157 + (qJD(1) * t84 + t80) * t185) * t125) * t126 + (t84 * t147 * t125 + (-t84 * t155 + (t84 * t159 + ((0.2e1 * t125 * t181 - t85 * t142 - t189 * t161) * t108 + (t118 * t126 * t145 + t125 * t85 + (t117 * t149 - (t85 - 0.2e1 * t162) * t125) * t111) * t109) * t179) * t125 + (-t91 + (-(t119 - t123) * t109 * t153 + t189 * t154) * t92) * t125 * qJD(1)) * t88) * t129, (-t164 * t184 + (-0.2e1 * t157 + (-qJD(2) * t83 - t80) * t185) * t129) * t128 + (t83 * t129 * t147 + (-t91 * t160 - ((-t126 * t81 - t163 * t98) * t109 + (t192 * t85 + t162 - t178) * t108) * t92 * t171 + (t129 * t159 + t92 * t164) * t83 - ((t81 - t163) * t108 + (t85 * t98 + qJD(2) + (-t85 - t178) * t126) * t109) * t128 * t179) * t88) * t125, 0, 0; 0.2e1 * (t101 * t138 + t105 * t176) * t186 + (0.2e1 * t105 * t156 - t144 * t101 * t168 + (t125 * t162 + t140) * t177 + (t138 * t87 + t105 * t86 - t140 * t175 - (qJD(2) * t125 * t127 + t144 * t124) * t106 * t126) * t102) * t96, t128 * t160 * t191 + (-t164 * t191 + (t139 * t158 + ((-qJD(3) * t101 - 0.2e1 * t156) * t127 + (-t127 * t86 + (-qJD(3) * t106 + t87) * t124) * t102) * t96) * t129) * t125, t158 + 0.2e1 * (-t86 * t102 * t96 + (-t102 * t186 - t96 * t183) * t106) * t106, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:19:24
	% EndTime: 2019-12-31 17:19:25
	% DurationCPUTime: 0.70s
	% Computational Cost: add. (1002->91), mult. (2519->210), div. (480->12), fcn. (2968->9), ass. (0->90)
	t131 = sin(qJ(1));
	t124 = t131 ^ 2;
	t130 = sin(qJ(2));
	t123 = t130 ^ 2;
	t133 = cos(qJ(2));
	t126 = 0.1e1 / t133 ^ 2;
	t180 = t123 * t126;
	t118 = t124 * t180 + 0.1e1;
	t122 = t130 * t123;
	t125 = 0.1e1 / t133;
	t179 = t125 * t130;
	t142 = qJD(2) * (t122 * t125 * t126 + t179);
	t134 = cos(qJ(1));
	t170 = qJD(1) * t134;
	t157 = t131 * t170;
	t186 = (t124 * t142 + t157 * t180) / t118 ^ 2;
	t196 = -0.2e1 * t186;
	t169 = qJD(2) * t131;
	t116 = 0.1e1 / t118;
	t156 = t126 * t169;
	t158 = t130 * t170;
	t168 = qJD(2) * t133;
	t90 = (-(-t131 * t168 - t158) * t125 + t123 * t156) * t116;
	t150 = t90 - t169;
	t153 = 0.1e1 + t180;
	t195 = t131 * t153;
	t129 = sin(qJ(3));
	t132 = cos(qJ(3));
	t172 = t133 * t134;
	t112 = t131 * t129 + t132 * t172;
	t151 = -t131 * t90 + qJD(2);
	t175 = t131 * t130;
	t115 = atan2(-t175, -t133);
	t114 = cos(t115);
	t113 = sin(t115);
	t161 = t113 * t175;
	t99 = -t114 * t133 - t161;
	t96 = 0.1e1 / t99;
	t106 = 0.1e1 / t112;
	t107 = 0.1e1 / t112 ^ 2;
	t97 = 0.1e1 / t99 ^ 2;
	t194 = t116 - 0.1e1;
	t128 = t134 ^ 2;
	t162 = t97 * t168;
	t181 = t114 * t130;
	t85 = t151 * t181 + (t150 * t133 - t158) * t113;
	t192 = t85 * t96 * t97;
	t95 = t123 * t128 * t97 + 0.1e1;
	t193 = (t128 * t130 * t162 + (-t128 * t192 - t97 * t157) * t123) / t95 ^ 2;
	t93 = 0.1e1 / t95;
	t191 = t93 * t97;
	t190 = t96 * t93;
	t174 = t131 * t132;
	t111 = t129 * t172 - t174;
	t105 = t111 ^ 2;
	t104 = t105 * t107 + 0.1e1;
	t183 = t107 * t111;
	t148 = -qJD(1) * t133 + qJD(3);
	t149 = qJD(3) * t133 - qJD(1);
	t167 = qJD(2) * t134;
	t155 = t130 * t167;
	t178 = t129 * t134;
	t92 = -t149 * t178 + (t148 * t131 - t155) * t132;
	t188 = t106 * t107 * t92;
	t173 = t131 * t133;
	t143 = t129 * t173 + t132 * t134;
	t91 = t143 * qJD(1) - qJD(3) * t112 + t129 * t155;
	t189 = 0.1e1 / t104 ^ 2 * (-t105 * t188 - t91 * t183);
	t184 = t106 * t129;
	t182 = t111 * t132;
	t177 = t130 * t134;
	t171 = qJD(1) * t131;
	t166 = 0.2e1 * t192;
	t165 = -0.2e1 * t189;
	t164 = t96 * t193;
	t163 = t111 * t188;
	t160 = t116 * t123 * t125;
	t154 = 0.2e1 * t97 * t193;
	t152 = t125 * t196;
	t147 = t131 * t160;
	t146 = t153 * t134;
	t145 = t148 * t134;
	t144 = t107 * t182 - t184;
	t110 = -t132 * t173 + t178;
	t103 = t116 * t195;
	t101 = 0.1e1 / t104;
	t89 = (t194 * t130 * t113 - t114 * t147) * t134;
	t88 = -t113 * t173 + t181 + (t113 * t133 - t114 * t175) * t103;
	t86 = t195 * t196 + (qJD(1) * t146 + 0.2e1 * t131 * t142) * t116;
	t1 = [t152 * t177 + (qJD(2) * t146 - t171 * t179) * t116, t86, 0, 0; (-t168 * t190 + (0.2e1 * t164 + (qJD(1) * t89 + t85) * t191) * t130) * t131 + (t89 * t154 * t130 + (-t89 * t162 + (t89 * t166 + ((0.2e1 * t130 * t186 - t90 * t147 - t194 * t168) * t113 + (t123 * t131 * t152 + t90 * t130 + (t122 * t156 - (t90 - 0.2e1 * t169) * t130) * t116) * t114) * t97 * t134) * t130 + (-t96 + (-(t124 - t128) * t114 * t160 + t194 * t161) * t97) * t130 * qJD(1)) * t93) * t134, (-t171 * t190 + (-0.2e1 * t164 + (-qJD(2) * t88 - t85) * t191) * t134) * t133 + (t88 * t134 * t154 + (-t96 * t167 - ((-t103 * t170 - t131 * t86) * t114 + (-t103 * t151 - t150) * t113) * t97 * t177 + (t134 * t166 + t97 * t171) * t88) * t93 - ((t86 - t170) * t113 + (t150 * t103 + t151) * t114) * t172 * t191) * t130, 0, 0; 0.2e1 * (t106 * t143 + t110 * t183) * t189 + (0.2e1 * t110 * t163 - t149 * t106 * t174 + (t130 * t169 + t145) * t184 + (t143 * t92 + t110 * t91 - t145 * t182 - (qJD(2) * t130 * t132 + t149 * t129) * t111 * t131) * t107) * t101, t144 * t165 * t177 + (t144 * t133 * t167 + (-t144 * t171 + ((-qJD(3) * t106 - 0.2e1 * t163) * t132 + (-t132 * t91 + (-qJD(3) * t111 + t92) * t129) * t107) * t134) * t130) * t101, t165 + 0.2e1 * (-t101 * t107 * t91 + (-t101 * t188 - t107 * t189) * t111) * t111, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end