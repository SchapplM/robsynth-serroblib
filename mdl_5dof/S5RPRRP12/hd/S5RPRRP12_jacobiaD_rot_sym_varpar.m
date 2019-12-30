% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRP12
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
%   Wie in S5RPRRP12_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 17:35
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRRP12_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP12_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_jacobiaD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:35:27
	% EndTime: 2019-12-29 17:35:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:35:27
	% EndTime: 2019-12-29 17:35:27
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:35:27
	% EndTime: 2019-12-29 17:35:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:35:27
	% EndTime: 2019-12-29 17:35:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:35:27
	% EndTime: 2019-12-29 17:35:28
	% DurationCPUTime: 1.38s
	% Computational Cost: add. (813->90), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->96)
	t132 = cos(qJ(1));
	t195 = 0.2e1 * t132;
	t167 = qJD(3) * t132;
	t126 = t132 ^ 2;
	t128 = sin(qJ(3));
	t121 = 0.1e1 / t128 ^ 2;
	t131 = cos(qJ(3));
	t125 = t131 ^ 2;
	t179 = t121 * t125;
	t118 = t126 * t179 + 0.1e1;
	t116 = 0.1e1 / t118;
	t120 = 0.1e1 / t128;
	t158 = t121 * t167;
	t129 = sin(qJ(1));
	t171 = qJD(1) * t131;
	t160 = t129 * t171;
	t90 = ((t128 * t167 + t160) * t120 + t125 * t158) * t116;
	t152 = -t90 + t167;
	t173 = t132 * t131;
	t115 = atan2(-t173, t128);
	t113 = sin(t115);
	t114 = cos(t115);
	t99 = -t113 * t173 + t114 * t128;
	t96 = 0.1e1 / t99;
	t127 = sin(qJ(4));
	t175 = t132 * t127;
	t130 = cos(qJ(4));
	t177 = t129 * t130;
	t110 = t128 * t177 + t175;
	t106 = 0.1e1 / t110;
	t107 = 0.1e1 / t110 ^ 2;
	t97 = 0.1e1 / t99 ^ 2;
	t194 = t116 - 0.1e1;
	t123 = t129 ^ 2;
	t170 = qJD(1) * t132;
	t147 = t125 * t129 * t170;
	t169 = qJD(3) * t128;
	t163 = t97 * t169;
	t153 = -t132 * t90 + qJD(3);
	t181 = t114 * t131;
	t85 = t153 * t181 + (t152 * t128 + t160) * t113;
	t192 = t85 * t96 * t97;
	t95 = t123 * t125 * t97 + 0.1e1;
	t193 = (t97 * t147 + (-t125 * t192 - t131 * t163) * t123) / t95 ^ 2;
	t93 = 0.1e1 / t95;
	t191 = t93 * t97;
	t190 = t96 * t93;
	t174 = t132 * t130;
	t178 = t129 * t127;
	t109 = t128 * t178 - t174;
	t105 = t109 ^ 2;
	t104 = t105 * t107 + 0.1e1;
	t184 = t107 * t109;
	t150 = qJD(1) * t128 + qJD(4);
	t143 = t150 * t132;
	t151 = qJD(4) * t128 + qJD(1);
	t145 = t151 * t127;
	t168 = qJD(3) * t131;
	t92 = t130 * t143 + (t130 * t168 - t145) * t129;
	t188 = t106 * t107 * t92;
	t144 = t151 * t130;
	t91 = t129 * t144 + (t129 * t168 + t143) * t127;
	t189 = 0.1e1 / t104 ^ 2 * (-t105 * t188 + t91 * t184);
	t187 = t129 * t97;
	t124 = t131 * t125;
	t141 = qJD(3) * (-t121 * t124 - t131) * t120;
	t186 = (-t121 * t147 + t126 * t141) / t118 ^ 2;
	t185 = t106 * t127;
	t183 = t109 * t130;
	t182 = t113 * t132;
	t180 = t120 * t125;
	t176 = t129 * t131;
	t172 = qJD(1) * t129;
	t166 = -0.2e1 * t192;
	t165 = 0.2e1 * t189;
	t164 = t96 * t193;
	t162 = t131 * t186;
	t161 = t116 * t180;
	t159 = t131 * t170;
	t157 = -0.2e1 * t97 * t193;
	t156 = 0.1e1 + t179;
	t155 = 0.2e1 * t109 * t188;
	t154 = t186 * t195;
	t149 = t132 * t161;
	t148 = t194 * t131 * t113;
	t146 = t156 * t129;
	t142 = t107 * t183 - t185;
	t140 = -t150 * t129 + t131 * t167;
	t112 = t128 * t174 - t178;
	t111 = t128 * t175 + t177;
	t103 = t156 * t132 * t116;
	t101 = 0.1e1 / t104;
	t89 = (-t114 * t149 - t148) * t129;
	t88 = t128 * t182 + t181 + (-t113 * t128 - t114 * t173) * t103;
	t86 = -t156 * t154 + (-qJD(1) * t146 + t141 * t195) * t116;
	t1 = [-0.2e1 * t129 * t120 * t162 + (-qJD(3) * t146 + t120 * t159) * t116, 0, t86, 0, 0; (t169 * t190 + (0.2e1 * t164 + (qJD(1) * t89 + t85) * t191) * t131) * t132 + (t89 * t157 * t131 + (-t89 * t163 + (t89 * t166 + ((t90 * t149 + t194 * t169 + 0.2e1 * t162) * t113 + (t154 * t180 + t131 * t90 + (t124 * t158 + (-t90 + 0.2e1 * t167) * t131) * t116) * t114) * t187) * t131 + (t96 + ((t123 - t126) * t114 * t161 - t132 * t148) * t97) * t171) * t93) * t129, 0, (t170 * t190 + (-0.2e1 * t164 + (-qJD(3) * t88 - t85) * t191) * t129) * t128 + (t88 * t129 * t157 + (t129 * qJD(3) * t96 + (-t114 * t132 * t86 + t152 * t113 + (-qJD(3) * t113 + t114 * t172 + t182 * t90) * t103) * t97 * t176 + (t129 * t166 + t97 * t170) * t88 + ((-t86 - t172) * t113 + (t152 * t103 - t153) * t114) * t128 * t187) * t93) * t131, 0, 0; (-t106 * t111 + t112 * t184) * t165 + (t112 * t155 + t132 * t106 * t144 + t140 * t185 + (t132 * t109 * t145 - t111 * t92 - t112 * t91 - t140 * t183) * t107) * t101, 0, t142 * t165 * t176 + (-t142 * t159 + (t142 * t169 + ((qJD(4) * t106 + t155) * t130 + (-t130 * t91 + (qJD(4) * t109 - t92) * t127) * t107) * t131) * t129) * t101, -0.2e1 * t189 + 0.2e1 * (t91 * t107 * t101 + (-t101 * t188 - t107 * t189) * t109) * t109, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:35:27
	% EndTime: 2019-12-29 17:35:28
	% DurationCPUTime: 1.38s
	% Computational Cost: add. (813->91), mult. (2519->208), div. (480->12), fcn. (2968->9), ass. (0->94)
	t135 = cos(qJ(1));
	t197 = 0.2e1 * t135;
	t131 = sin(qJ(3));
	t134 = cos(qJ(3));
	t176 = t135 * t134;
	t118 = atan2(-t176, t131);
	t116 = sin(t118);
	t117 = cos(t118);
	t102 = -t116 * t176 + t117 * t131;
	t99 = 0.1e1 / t102;
	t130 = sin(qJ(4));
	t132 = sin(qJ(1));
	t133 = cos(qJ(4));
	t179 = t132 * t133;
	t113 = t130 * t135 + t131 * t179;
	t109 = 0.1e1 / t113;
	t123 = 0.1e1 / t131;
	t100 = 0.1e1 / t102 ^ 2;
	t110 = 0.1e1 / t113 ^ 2;
	t124 = 0.1e1 / t131 ^ 2;
	t129 = t135 ^ 2;
	t128 = t134 ^ 2;
	t183 = t124 * t128;
	t121 = t129 * t183 + 0.1e1;
	t119 = 0.1e1 / t121;
	t196 = t119 - 0.1e1;
	t126 = t132 ^ 2;
	t173 = qJD(1) * t135;
	t150 = t128 * t132 * t173;
	t171 = qJD(3) * t134;
	t182 = t126 * t128;
	t170 = qJD(3) * t135;
	t160 = t124 * t170;
	t174 = qJD(1) * t134;
	t162 = t132 * t174;
	t93 = ((t131 * t170 + t162) * t123 + t128 * t160) * t119;
	t155 = -t93 + t170;
	t156 = -t135 * t93 + qJD(3);
	t185 = t117 * t134;
	t88 = t156 * t185 + (t155 * t131 + t162) * t116;
	t193 = t99 * t100 * t88;
	t98 = t100 * t182 + 0.1e1;
	t195 = (-t182 * t193 + (-t126 * t131 * t171 + t150) * t100) / t98 ^ 2;
	t96 = 0.1e1 / t98;
	t194 = t100 * t96;
	t177 = t135 * t133;
	t180 = t132 * t130;
	t112 = t131 * t180 - t177;
	t108 = t112 ^ 2;
	t107 = t108 * t110 + 0.1e1;
	t187 = t110 * t112;
	t153 = qJD(1) * t131 + qJD(4);
	t146 = t153 * t135;
	t154 = qJD(4) * t131 + qJD(1);
	t148 = t154 * t130;
	t95 = t133 * t146 + (t133 * t171 - t148) * t132;
	t191 = t109 * t110 * t95;
	t147 = t154 * t133;
	t94 = t132 * t147 + (t132 * t171 + t146) * t130;
	t192 = 0.1e1 / t107 ^ 2 * (-t108 * t191 + t94 * t187);
	t127 = t134 * t128;
	t144 = qJD(3) * (-t124 * t127 - t134) * t123;
	t189 = (-t124 * t150 + t129 * t144) / t121 ^ 2;
	t188 = t109 * t130;
	t186 = t112 * t133;
	t184 = t123 * t128;
	t181 = t131 * t135;
	t178 = t132 * t134;
	t175 = qJD(1) * t132;
	t172 = qJD(3) * t131;
	t169 = -0.2e1 * t195;
	t168 = -0.2e1 * t193;
	t167 = 0.2e1 * t192;
	t166 = t99 * t195;
	t165 = t96 * t172;
	t164 = t134 * t189;
	t163 = t119 * t184;
	t161 = t134 * t173;
	t159 = 0.1e1 + t183;
	t158 = 0.2e1 * t112 * t191;
	t157 = t189 * t197;
	t152 = t135 * t163;
	t151 = t196 * t134 * t116;
	t149 = t159 * t132;
	t145 = t110 * t186 - t188;
	t143 = -t153 * t132 + t134 * t170;
	t115 = t131 * t177 - t180;
	t114 = t130 * t181 + t179;
	t106 = t159 * t135 * t119;
	t104 = 0.1e1 / t107;
	t92 = (-t117 * t152 - t151) * t132;
	t91 = t116 * t181 + t185 + (-t116 * t131 - t117 * t176) * t106;
	t89 = -t159 * t157 + (-qJD(1) * t149 + t144 * t197) * t119;
	t1 = [-0.2e1 * t132 * t123 * t164 + (-qJD(3) * t149 + t123 * t161) * t119, 0, t89, 0, 0; (t99 * t165 + (0.2e1 * t166 + (qJD(1) * t92 + t88) * t194) * t134) * t135 + (t92 * t134 * t96 * t168 + (-t92 * t165 + (t92 * t169 + ((t93 * t152 + t196 * t172 + 0.2e1 * t164) * t116 + (t157 * t184 + t134 * t93 + (t127 * t160 + (-t93 + 0.2e1 * t170) * t134) * t119) * t117) * t96 * t132) * t134) * t100 + (t99 + ((t126 - t129) * t117 * t163 - t135 * t151) * t100) * t96 * t174) * t132, 0, (t99 * t96 * t173 + (-0.2e1 * t166 + (-qJD(3) * t91 - t88) * t194) * t132) * t131 + (((qJD(3) * t99 + t91 * t168) * t132 + (t91 * t173 + ((t106 * t175 - t135 * t89) * t117 + ((t106 * t135 - 0.1e1) * t93 + (-t106 + t135) * qJD(3)) * t116) * t178) * t100) * t96 + (t91 * t169 + ((-t89 - t175) * t116 + (t155 * t106 - t156) * t117) * t96 * t131) * t100 * t132) * t134, 0, 0; (-t109 * t114 + t115 * t187) * t167 + (t115 * t158 + t135 * t109 * t147 + t143 * t188 + (t135 * t112 * t148 - t114 * t95 - t115 * t94 - t143 * t186) * t110) * t104, 0, t145 * t167 * t178 + (-t145 * t161 + (t145 * t172 + ((qJD(4) * t109 + t158) * t133 + (-t133 * t94 + (qJD(4) * t112 - t95) * t130) * t110) * t134) * t132) * t104, -0.2e1 * t192 + 0.2e1 * (t104 * t110 * t94 + (-t104 * t191 - t110 * t192) * t112) * t112, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end