% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRP5
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
%   Wie in S6RPPRRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:54
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRP5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:54:38
	% EndTime: 2019-10-09 23:54:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:54:38
	% EndTime: 2019-10-09 23:54:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:54:38
	% EndTime: 2019-10-09 23:54:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:54:38
	% EndTime: 2019-10-09 23:54:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:54:38
	% EndTime: 2019-10-09 23:54:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:54:38
	% EndTime: 2019-10-09 23:54:39
	% DurationCPUTime: 0.95s
	% Computational Cost: add. (624->91), mult. (2519->209), div. (480->12), fcn. (2968->9), ass. (0->95)
	t131 = sin(qJ(1));
	t171 = qJD(4) * t131;
	t125 = t131 ^ 2;
	t130 = sin(qJ(4));
	t123 = 0.1e1 / t130 ^ 2;
	t133 = cos(qJ(4));
	t127 = t133 ^ 2;
	t183 = t123 * t127;
	t118 = t125 * t183 + 0.1e1;
	t115 = 0.1e1 / t118;
	t122 = 0.1e1 / t130;
	t158 = t123 * t171;
	t134 = cos(qJ(1));
	t173 = qJD(1) * t134;
	t159 = t133 * t173;
	t90 = ((-t130 * t171 + t159) * t122 - t127 * t158) * t115;
	t150 = -t90 - t171;
	t154 = 0.1e1 + t183;
	t196 = t131 * t154;
	t178 = t131 * t133;
	t117 = atan2(t178, t130);
	t114 = cos(t117);
	t113 = sin(t117);
	t163 = t113 * t178;
	t99 = t114 * t130 + t163;
	t96 = 0.1e1 / t99;
	t132 = cos(qJ(5));
	t177 = t132 * t134;
	t161 = t130 * t177;
	t129 = sin(qJ(5));
	t180 = t131 * t129;
	t112 = t161 - t180;
	t106 = 0.1e1 / t112;
	t107 = 0.1e1 / t112 ^ 2;
	t97 = 0.1e1 / t99 ^ 2;
	t195 = t115 - 0.1e1;
	t151 = t131 * t90 + qJD(4);
	t184 = t114 * t133;
	t85 = t151 * t184 + (t150 * t130 + t159) * t113;
	t194 = t85 * t96 * t97;
	t128 = t134 ^ 2;
	t95 = t127 * t128 * t97 + 0.1e1;
	t93 = 0.1e1 / t95;
	t193 = t93 * t97;
	t192 = t96 * t93;
	t179 = t131 * t132;
	t181 = t130 * t134;
	t111 = t129 * t181 + t179;
	t105 = t111 ^ 2;
	t104 = t105 * t107 + 0.1e1;
	t187 = t107 * t111;
	t148 = qJD(1) * t130 + qJD(5);
	t149 = qJD(5) * t130 + qJD(1);
	t169 = qJD(4) * t134;
	t157 = t133 * t169;
	t182 = t129 * t134;
	t92 = -t149 * t182 + (-t148 * t131 + t157) * t132;
	t190 = t106 * t107 * t92;
	t91 = -qJD(5) * t161 - t129 * t157 - t132 * t173 + t148 * t180;
	t191 = 0.1e1 / t104 ^ 2 * (-t105 * t190 - t91 * t187);
	t126 = t133 * t127;
	t145 = (t123 * t126 + t133) * t122;
	t160 = t131 * t173;
	t189 = (-t145 * t125 * qJD(4) + t160 * t183) / t118 ^ 2;
	t188 = t106 * t129;
	t186 = t111 * t132;
	t185 = t113 * t131;
	t176 = t133 * t134;
	t175 = qJD(1) * t131;
	t174 = qJD(1) * t133;
	t172 = qJD(4) * t130;
	t170 = qJD(4) * t133;
	t164 = t97 * t172;
	t168 = -0.2e1 * (-t128 * t133 * t164 + (-t128 * t194 - t97 * t160) * t127) / t95 ^ 2;
	t167 = -0.2e1 * t194;
	t166 = 0.2e1 * t191;
	t165 = 0.2e1 * t189;
	t162 = t115 * t122 * t127;
	t156 = t96 * t168;
	t155 = t97 * t168;
	t153 = 0.2e1 * t111 * t190;
	t152 = -0.2e1 * t122 * t189;
	t147 = t131 * t162;
	t146 = t154 * t134;
	t144 = t148 * t134;
	t143 = t107 * t186 - t188;
	t142 = t143 * t134;
	t110 = -t130 * t179 - t182;
	t109 = -t130 * t180 + t177;
	t103 = t115 * t196;
	t101 = 0.1e1 / t104;
	t89 = (-t195 * t133 * t113 + t114 * t147) * t134;
	t88 = -t130 * t185 + t184 - (-t113 * t130 + t114 * t178) * t103;
	t86 = t165 * t196 + (-qJD(1) * t146 + 0.2e1 * t145 * t171) * t115;
	t1 = [t152 * t176 + (-t122 * t131 * t174 - qJD(4) * t146) * t115, 0, 0, t86, 0, 0; (-t172 * t192 + (t156 + (-qJD(1) * t89 - t85) * t193) * t133) * t131 + (t89 * t155 * t133 + (-t89 * t164 + (t89 * t167 + ((t133 * t165 - t90 * t147 + t195 * t172) * t113 + (t127 * t131 * t152 + t133 * t90 + (-t126 * t158 + (-t90 - 0.2e1 * t171) * t133) * t115) * t114) * t97 * t134) * t133 + (t96 + ((-t125 + t128) * t114 * t162 + t195 * t163) * t97) * t174) * t93) * t134, 0, 0, (-t175 * t192 + (t156 + (-qJD(4) * t88 - t85) * t193) * t134) * t130 + (t88 * t134 * t155 + (t96 * t169 + (t114 * t131 * t86 + t150 * t113 - (-qJD(4) * t113 + t114 * t173 - t185 * t90) * t103) * t97 * t176 + (t134 * t167 - t97 * t175) * t88) * t93 + ((-t86 - t173) * t113 + (-t150 * t103 - t151) * t114) * t181 * t193) * t133, 0, 0; (-t106 * t109 + t110 * t187) * t166 + (t110 * t153 - t149 * t106 * t179 + (-t131 * t170 - t144) * t188 + (-t109 * t92 + t110 * t91 + t144 * t186 - (t149 * t129 - t132 * t170) * t111 * t131) * t107) * t101, 0, 0, t133 * t142 * t166 + (t142 * t172 + (t143 * t175 + ((qJD(5) * t106 + t153) * t132 + (t132 * t91 + (qJD(5) * t111 - t92) * t129) * t107) * t134) * t133) * t101, -0.2e1 * t191 + 0.2e1 * (-t101 * t107 * t91 + (-t101 * t190 - t107 * t191) * t111) * t111, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:54:38
	% EndTime: 2019-10-09 23:54:39
	% DurationCPUTime: 0.94s
	% Computational Cost: add. (624->92), mult. (2519->208), div. (480->12), fcn. (2968->9), ass. (0->94)
	t134 = sin(qJ(1));
	t133 = sin(qJ(4));
	t126 = 0.1e1 / t133 ^ 2;
	t136 = cos(qJ(4));
	t130 = t136 ^ 2;
	t186 = t126 * t130;
	t158 = 0.1e1 + t186;
	t198 = t134 * t158;
	t179 = t134 * t136;
	t120 = atan2(t179, t133);
	t117 = cos(t120);
	t116 = sin(t120);
	t165 = t116 * t179;
	t102 = t117 * t133 + t165;
	t99 = 0.1e1 / t102;
	t135 = cos(qJ(5));
	t137 = cos(qJ(1));
	t178 = t135 * t137;
	t163 = t133 * t178;
	t132 = sin(qJ(5));
	t181 = t134 * t132;
	t115 = t163 - t181;
	t109 = 0.1e1 / t115;
	t125 = 0.1e1 / t133;
	t100 = 0.1e1 / t102 ^ 2;
	t110 = 0.1e1 / t115 ^ 2;
	t128 = t134 ^ 2;
	t121 = t128 * t186 + 0.1e1;
	t118 = 0.1e1 / t121;
	t197 = t118 - 0.1e1;
	t131 = t137 ^ 2;
	t185 = t130 * t131;
	t98 = t100 * t185 + 0.1e1;
	t96 = 0.1e1 / t98;
	t196 = t100 * t96;
	t172 = qJD(4) * t134;
	t161 = t126 * t172;
	t174 = qJD(1) * t137;
	t162 = t136 * t174;
	t93 = ((-t133 * t172 + t162) * t125 - t130 * t161) * t118;
	t154 = -t93 - t172;
	t155 = t134 * t93 + qJD(4);
	t187 = t117 * t136;
	t88 = t155 * t187 + (t154 * t133 + t162) * t116;
	t195 = t99 * t100 * t88;
	t180 = t134 * t135;
	t182 = t133 * t137;
	t114 = t132 * t182 + t180;
	t108 = t114 ^ 2;
	t107 = t108 * t110 + 0.1e1;
	t190 = t110 * t114;
	t152 = qJD(1) * t133 + qJD(5);
	t153 = qJD(5) * t133 + qJD(1);
	t171 = qJD(4) * t136;
	t160 = t137 * t171;
	t183 = t132 * t137;
	t95 = -t153 * t183 + (-t152 * t134 + t160) * t135;
	t193 = t109 * t110 * t95;
	t94 = -qJD(5) * t163 - t132 * t160 - t135 * t174 + t152 * t181;
	t194 = 0.1e1 / t107 ^ 2 * (-t108 * t193 - t94 * t190);
	t129 = t136 * t130;
	t148 = (t126 * t129 + t136) * t125;
	t184 = t130 * t134;
	t150 = t174 * t184;
	t192 = (-t148 * t128 * qJD(4) + t126 * t150) / t121 ^ 2;
	t191 = t109 * t132;
	t189 = t114 * t135;
	t188 = t116 * t133;
	t177 = t136 * t137;
	t176 = qJD(1) * t134;
	t175 = qJD(1) * t136;
	t173 = qJD(4) * t133;
	t170 = -0.2e1 * (-t185 * t195 + (-t131 * t133 * t171 - t150) * t100) / t98 ^ 2;
	t169 = -0.2e1 * t195;
	t168 = 0.2e1 * t194;
	t167 = 0.2e1 * t192;
	t166 = t96 * t173;
	t164 = t118 * t125 * t130;
	t159 = t99 * t170;
	t157 = 0.2e1 * t114 * t193;
	t156 = -0.2e1 * t125 * t192;
	t151 = t134 * t164;
	t149 = t158 * t137;
	t147 = t152 * t137;
	t146 = t110 * t189 - t191;
	t145 = t146 * t137;
	t113 = -t133 * t180 - t183;
	t112 = -t133 * t181 + t178;
	t106 = t118 * t198;
	t104 = 0.1e1 / t107;
	t92 = (-t197 * t136 * t116 + t117 * t151) * t137;
	t91 = -t134 * t188 + t187 - (t117 * t179 - t188) * t106;
	t89 = t167 * t198 + (-qJD(1) * t149 + 0.2e1 * t148 * t172) * t118;
	t1 = [t156 * t177 + (-t125 * t134 * t175 - qJD(4) * t149) * t118, 0, 0, t89, 0, 0; (-t99 * t166 + (t159 + (-qJD(1) * t92 - t88) * t196) * t136) * t134 + (t92 * t136 * t96 * t169 + (-t92 * t166 + (t92 * t170 + ((t136 * t167 - t93 * t151 + t197 * t173) * t116 + (t156 * t184 + t136 * t93 + (-t129 * t161 + (-t93 - 0.2e1 * t172) * t136) * t118) * t117) * t96 * t137) * t136) * t100 + (t99 + ((-t128 + t131) * t117 * t164 + t197 * t165) * t100) * t96 * t175) * t137, 0, 0, (-t99 * t96 * t176 + (t159 + (-qJD(4) * t91 - t88) * t196) * t137) * t133 + (t91 * t137 * t100 * t170 + ((qJD(4) * t99 + t91 * t169) * t137 + (-t91 * t176 + ((-t106 * t174 + t134 * t89) * t117 + ((t106 * t134 - 0.1e1) * t93 + (t106 - t134) * qJD(4)) * t116) * t177) * t100) * t96 + ((-t89 - t174) * t116 + (-t154 * t106 - t155) * t117) * t182 * t196) * t136, 0, 0; (-t109 * t112 + t113 * t190) * t168 + (t113 * t157 - t153 * t109 * t180 + (-t134 * t171 - t147) * t191 + (-t112 * t95 + t113 * t94 + t147 * t189 - (t153 * t132 - t135 * t171) * t114 * t134) * t110) * t104, 0, 0, t136 * t145 * t168 + (t145 * t173 + (t146 * t176 + ((qJD(5) * t109 + t157) * t135 + (t135 * t94 + (qJD(5) * t114 - t95) * t132) * t110) * t137) * t136) * t104, -0.2e1 * t194 + 0.2e1 * (-t104 * t110 * t94 + (-t104 * t193 - t110 * t194) * t114) * t114, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end