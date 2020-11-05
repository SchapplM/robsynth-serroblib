% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR13 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:12
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRPR13_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR13_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:24
	% EndTime: 2020-11-04 22:12:24
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:24
	% EndTime: 2020-11-04 22:12:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t91 = cos(qJ(1));
	t90 = sin(qJ(1));
	t1 = [t91, -t90, 0, 0; t90, t91, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:24
	% EndTime: 2020-11-04 22:12:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t92 = sin(pkin(6));
	t95 = sin(qJ(1));
	t103 = t95 * t92;
	t94 = sin(qJ(2));
	t102 = t95 * t94;
	t96 = cos(qJ(2));
	t101 = t95 * t96;
	t97 = cos(qJ(1));
	t100 = t97 * t92;
	t99 = t97 * t94;
	t98 = t97 * t96;
	t93 = cos(pkin(6));
	t1 = [-t93 * t102 + t98, -t93 * t101 - t99, t103, t97 * pkin(1) + pkin(8) * t103 + 0; t93 * t99 + t101, t93 * t98 - t102, -t100, t95 * pkin(1) - pkin(8) * t100 + 0; t92 * t94, t92 * t96, t93, t93 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:24
	% EndTime: 2020-11-04 22:12:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (23->18), mult. (39->24), div. (0->0), fcn. (52->6), ass. (0->14)
	t108 = sin(qJ(2));
	t109 = sin(qJ(1));
	t116 = t109 * t108;
	t110 = cos(qJ(2));
	t115 = t109 * t110;
	t111 = cos(qJ(1));
	t114 = t111 * t108;
	t113 = t111 * t110;
	t112 = pkin(2) * t108 - qJ(3) * t110;
	t107 = cos(pkin(6));
	t106 = sin(pkin(6));
	t105 = t110 * pkin(2) + t108 * qJ(3) + pkin(1);
	t104 = t106 * pkin(8) - t112 * t107;
	t1 = [t109 * t106, t107 * t116 - t113, t107 * t115 + t114, t104 * t109 + t105 * t111 + 0; -t111 * t106, -t107 * t114 - t115, -t107 * t113 + t116, -t104 * t111 + t105 * t109 + 0; t107, -t106 * t108, -t106 * t110, t107 * pkin(8) + t112 * t106 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:24
	% EndTime: 2020-11-04 22:12:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (36->24), mult. (60->37), div. (0->0), fcn. (81->8), ass. (0->22)
	t120 = sin(pkin(6));
	t125 = cos(qJ(4));
	t137 = t120 * t125;
	t122 = sin(qJ(4));
	t136 = t122 * t120;
	t123 = sin(qJ(2));
	t124 = sin(qJ(1));
	t135 = t124 * t123;
	t126 = cos(qJ(2));
	t134 = t124 * t126;
	t133 = t125 * t126;
	t127 = cos(qJ(1));
	t132 = t127 * t123;
	t131 = t127 * t126;
	t129 = pkin(2) + pkin(9);
	t130 = qJ(3) * t126 - t123 * t129;
	t128 = pkin(3) + pkin(8);
	t121 = cos(pkin(6));
	t119 = t123 * qJ(3) + t129 * t126 + pkin(1);
	t118 = -t121 * t131 + t135;
	t117 = t120 * t128 + t130 * t121;
	t1 = [t124 * t137 + (t121 * t134 + t132) * t122, (t121 * t133 - t136) * t124 + t125 * t132, -t121 * t135 + t131, t117 * t124 + t119 * t127 + 0; t118 * t122 - t127 * t137, t118 * t125 + t127 * t136, t121 * t132 + t134, -t117 * t127 + t119 * t124 + 0; t121 * t125 - t126 * t136, -t120 * t133 - t121 * t122, t120 * t123, -t130 * t120 + t128 * t121 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:24
	% EndTime: 2020-11-04 22:12:24
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (65->34), mult. (116->57), div. (0->0), fcn. (150->10), ass. (0->29)
	t147 = sin(pkin(6));
	t150 = sin(qJ(4));
	t165 = t147 * t150;
	t151 = sin(qJ(2));
	t164 = t147 * t151;
	t149 = cos(pkin(6));
	t163 = t149 * t151;
	t162 = t150 * t151;
	t154 = cos(qJ(2));
	t161 = t150 * t154;
	t153 = cos(qJ(4));
	t160 = t153 * t154;
	t145 = pkin(4) * t150 - qJ(5) * t153 + qJ(3);
	t156 = pkin(2) + pkin(9);
	t159 = t145 * t154 - t151 * t156;
	t158 = t147 * t153 + t149 * t161;
	t157 = pkin(4) * t153 + qJ(5) * t150 + pkin(3) + pkin(8);
	t155 = cos(qJ(1));
	t152 = sin(qJ(1));
	t148 = cos(pkin(11));
	t146 = sin(pkin(11));
	t144 = -t147 * t161 + t149 * t153;
	t143 = t146 * t154 + t148 * t162;
	t142 = t146 * t162 - t148 * t154;
	t141 = t145 * t151 + t154 * t156 + pkin(1);
	t140 = -t146 * t163 + t158 * t148;
	t139 = t158 * t146 + t148 * t163;
	t138 = t147 * t157 + t159 * t149;
	t1 = [t140 * t152 + t143 * t155, -t139 * t152 - t142 * t155, (-t149 * t160 + t165) * t152 - t153 * t151 * t155, t138 * t152 + t141 * t155 + 0; -t140 * t155 + t143 * t152, t139 * t155 - t142 * t152, (t149 * t154 * t155 - t151 * t152) * t153 - t155 * t165, -t138 * t155 + t141 * t152 + 0; t144 * t148 + t146 * t164, -t144 * t146 + t148 * t164, t147 * t160 + t149 * t150, -t159 * t147 + t157 * t149 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:12:24
	% EndTime: 2020-11-04 22:12:24
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (98->38), mult. (125->56), div. (0->0), fcn. (159->12), ass. (0->33)
	t178 = sin(pkin(6));
	t182 = sin(qJ(2));
	t197 = t178 * t182;
	t181 = sin(qJ(4));
	t196 = t181 * t178;
	t185 = cos(qJ(2));
	t195 = t181 * t185;
	t186 = cos(qJ(1));
	t194 = t182 * t186;
	t183 = sin(qJ(1));
	t193 = t183 * t182;
	t184 = cos(qJ(4));
	t192 = t184 * t185;
	t191 = t186 * t185;
	t174 = cos(pkin(11)) * pkin(5) + pkin(4);
	t180 = qJ(5) + pkin(10);
	t168 = t174 * t181 - t180 * t184 + qJ(3);
	t173 = sin(pkin(11)) * pkin(5) + pkin(2) + pkin(9);
	t190 = t168 * t185 - t173 * t182;
	t179 = cos(pkin(6));
	t169 = t184 * t178 + t179 * t195;
	t189 = t169 * t183 + t181 * t194;
	t188 = -t169 * t186 + t181 * t193;
	t187 = t174 * t184 + t180 * t181 + pkin(3) + pkin(8);
	t177 = pkin(11) + qJ(6);
	t176 = cos(t177);
	t175 = sin(t177);
	t172 = -t179 * t193 + t191;
	t171 = t179 * t194 + t183 * t185;
	t170 = -t178 * t195 + t179 * t184;
	t167 = t168 * t182 + t173 * t185 + pkin(1);
	t166 = t178 * t187 + t190 * t179;
	t1 = [t172 * t175 + t189 * t176, t172 * t176 - t189 * t175, (-t179 * t192 + t196) * t183 - t184 * t194, t166 * t183 + t167 * t186 + 0; t171 * t175 + t188 * t176, t171 * t176 - t188 * t175, (t179 * t191 - t193) * t184 - t186 * t196, -t166 * t186 + t167 * t183 + 0; t170 * t176 + t175 * t197, -t170 * t175 + t176 * t197, t178 * t192 + t179 * t181, -t190 * t178 + t187 * t179 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end