% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRR6
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
%   Wie in S5RRPRR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRPRR6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
JaD_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (74->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:06
	% EndTime: 2022-01-20 11:18:06
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (1358->41), mult. (1512->108), div. (234->12), fcn. (1884->9), ass. (0->57)
	t127 = qJ(1) + qJ(2);
	t120 = cos(t127);
	t118 = t120 ^ 2;
	t128 = sin(pkin(9));
	t122 = t128 ^ 2;
	t158 = t118 * t122;
	t119 = sin(t127);
	t117 = t119 ^ 2;
	t129 = cos(pkin(9));
	t124 = 0.1e1 / t129 ^ 2;
	t113 = t117 * t122 * t124 + 0.1e1;
	t111 = 0.1e1 / t113;
	t157 = (t111 - 0.1e1) * t128;
	t130 = sin(qJ(4));
	t131 = cos(qJ(4));
	t144 = t129 * t131;
	t107 = t119 * t130 + t120 * t144;
	t148 = t119 * t128;
	t110 = atan2(-t148, -t129);
	t108 = sin(t110);
	t109 = cos(t110);
	t96 = -t108 * t148 - t109 * t129;
	t93 = 0.1e1 / t96;
	t101 = 0.1e1 / t107;
	t123 = 0.1e1 / t129;
	t102 = 0.1e1 / t107 ^ 2;
	t94 = 0.1e1 / t96 ^ 2;
	t145 = t129 * t130;
	t106 = -t119 * t131 + t120 * t145;
	t100 = t106 ^ 2;
	t126 = qJD(1) + qJD(2);
	t140 = t119 * t145 + t120 * t131;
	t91 = -t107 * qJD(4) + t140 * t126;
	t151 = t91 * t102;
	t105 = -t119 * t144 + t120 * t130;
	t92 = -t106 * qJD(4) + t105 * t126;
	t154 = t101 * t102 * t92;
	t99 = t100 * t102 + 0.1e1;
	t155 = (-t100 * t154 - t106 * t151) / t99 ^ 2;
	t153 = t119 * t94;
	t152 = t120 * t94;
	t150 = t105 * t106;
	t149 = t119 * t126;
	t146 = t122 * t123;
	t143 = t111 * t146;
	t112 = 0.1e1 / t113 ^ 2;
	t142 = t112 * t128 * t158;
	t86 = (-t109 * t119 * t143 + t108 * t157) * t120;
	t125 = t123 * t124;
	t97 = 0.1e1 / t99;
	t95 = t93 * t94;
	t90 = t94 * t158 + 0.1e1;
	t87 = (-t111 * t123 * t128 - 0.2e1 * t125 * t142) * t149;
	t85 = t126 * t86;
	t82 = 0.2e1 * (t101 * t140 + t102 * t150) * t155 + ((t105 * qJD(4) - t106 * t126) * t101 + 0.2e1 * t150 * t154 + (t140 * t92 - (t140 * qJD(4) - t107 * t126) * t106 + t105 * t91) * t102) * t97;
	t81 = (0.2e1 * (t119 * t93 + t86 * t152) / t90 ^ 2 * (-t118 * t85 * t95 - t149 * t152) * t122 + ((0.2e1 * t120 * t86 * t95 + t153) * t85 + (t86 * t153 + (-t93 + (t124 * t142 + t157) * t108 * t153 - (t117 * t143 + (-0.2e1 * t143 + (0.2e1 * t117 * t122 ^ 2 * t125 + t146) * t112) * t118) * t94 * t109) * t120) * t126) / t90) * t128;
	t1 = [t87, t87, 0, 0, 0; t81, t81, 0, 0, 0; t82, t82, 0, -0.2e1 * t155 + 0.2e1 * (-t97 * t151 + (-t102 * t155 - t97 * t154) * t106) * t106, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:06
	% EndTime: 2022-01-20 11:18:06
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (1968->45), mult. (1706->112), div. (252->12), fcn. (2094->9), ass. (0->65)
	t158 = qJ(1) + qJ(2);
	t149 = cos(t158);
	t145 = t149 ^ 2;
	t159 = sin(pkin(9));
	t151 = t159 ^ 2;
	t191 = t145 * t151;
	t147 = sin(t158);
	t144 = t147 ^ 2;
	t160 = cos(pkin(9));
	t153 = 0.1e1 / t160 ^ 2;
	t142 = t144 * t151 * t153 + 0.1e1;
	t140 = 0.1e1 / t142;
	t190 = (t140 - 0.1e1) * t159;
	t157 = qJ(4) + qJ(5);
	t148 = cos(t157);
	t177 = t149 * t160;
	t146 = sin(t157);
	t183 = t147 * t146;
	t134 = t148 * t177 + t183;
	t180 = t147 * t159;
	t138 = atan2(-t180, -t160);
	t136 = sin(t138);
	t137 = cos(t138);
	t126 = -t136 * t180 - t137 * t160;
	t123 = 0.1e1 / t126;
	t128 = 0.1e1 / t134;
	t152 = 0.1e1 / t160;
	t124 = 0.1e1 / t126 ^ 2;
	t129 = 0.1e1 / t134 ^ 2;
	t182 = t147 * t148;
	t133 = t146 * t177 - t182;
	t127 = t133 ^ 2;
	t122 = t127 * t129 + 0.1e1;
	t155 = qJD(4) + qJD(5);
	t156 = qJD(1) + qJD(2);
	t171 = t155 * t160 - t156;
	t172 = -t156 * t160 + t155;
	t178 = t149 * t146;
	t115 = -t171 * t178 + t172 * t182;
	t186 = t115 * t128 * t129;
	t179 = t147 * t160;
	t169 = t146 * t179 + t149 * t148;
	t114 = -t134 * t155 + t169 * t156;
	t187 = t114 * t129;
	t188 = (-t127 * t186 - t133 * t187) / t122 ^ 2;
	t185 = t124 * t149;
	t132 = -t148 * t179 + t178;
	t184 = t132 * t133;
	t181 = t147 * t156;
	t176 = t151 * t152;
	t175 = t140 * t176;
	t141 = 0.1e1 / t142 ^ 2;
	t174 = t141 * t159 * t191;
	t170 = t172 * t149;
	t113 = (-t137 * t147 * t175 + t136 * t190) * t149;
	t154 = t152 * t153;
	t125 = t123 * t124;
	t120 = 0.1e1 / t122;
	t119 = t124 * t191 + 0.1e1;
	t116 = (-t140 * t152 * t159 - 0.2e1 * t154 * t174) * t181;
	t112 = t156 * t113;
	t109 = -0.2e1 * t188 + 0.2e1 * (-t120 * t187 + (-t120 * t186 - t129 * t188) * t133) * t133;
	t108 = 0.2e1 * (t128 * t169 + t129 * t184) * t188 + ((t146 * t170 - t171 * t182) * t128 + 0.2e1 * t184 * t186 + (t169 * t115 - (t148 * t170 + t171 * t183) * t133 + t132 * t114) * t129) * t120;
	t107 = (0.2e1 * (t113 * t185 + t123 * t147) / t119 ^ 2 * (-t112 * t125 * t145 - t181 * t185) * t151 + ((0.2e1 * t113 * t125 * t149 + t124 * t147) * t112 + (-t149 * t123 + ((t113 + (t153 * t174 + t190) * t149 * t136) * t147 - (t144 * t175 + (-0.2e1 * t175 + (0.2e1 * t144 * t151 ^ 2 * t154 + t176) * t141) * t145) * t149 * t137) * t124) * t156) / t119) * t159;
	t1 = [t116, t116, 0, 0, 0; t107, t107, 0, 0, 0; t108, t108, 0, t109, t109;];
	JaD_rot = t1;
end