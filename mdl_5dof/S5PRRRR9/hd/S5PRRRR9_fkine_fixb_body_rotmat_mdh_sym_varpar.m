% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR9 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:09
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:43
	% EndTime: 2020-11-04 20:09:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:43
	% EndTime: 2020-11-04 20:09:43
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t88 = cos(pkin(10));
	t87 = sin(pkin(10));
	t1 = [t88, -t87, 0, 0; t87, t88, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:43
	% EndTime: 2020-11-04 20:09:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t89 = sin(pkin(10));
	t90 = sin(pkin(5));
	t98 = t89 * t90;
	t91 = cos(pkin(10));
	t97 = t91 * t90;
	t92 = cos(pkin(5));
	t93 = sin(qJ(2));
	t96 = t92 * t93;
	t94 = cos(qJ(2));
	t95 = t92 * t94;
	t1 = [-t89 * t96 + t91 * t94, -t89 * t95 - t91 * t93, t98, t91 * pkin(1) + pkin(6) * t98 + 0; t89 * t94 + t91 * t96, -t89 * t93 + t91 * t95, -t97, t89 * pkin(1) - pkin(6) * t97 + 0; t90 * t93, t90 * t94, t92, t92 * pkin(6) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:43
	% EndTime: 2020-11-04 20:09:43
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t102 = sin(pkin(5));
	t115 = t102 * pkin(6);
	t101 = sin(pkin(10));
	t104 = cos(pkin(5));
	t114 = t101 * t104;
	t105 = sin(qJ(3));
	t113 = t102 * t105;
	t107 = cos(qJ(3));
	t112 = t102 * t107;
	t103 = cos(pkin(10));
	t111 = t103 * t104;
	t106 = sin(qJ(2));
	t110 = t104 * t106;
	t108 = cos(qJ(2));
	t109 = t104 * t108;
	t100 = -t101 * t110 + t103 * t108;
	t99 = t101 * t108 + t103 * t110;
	t1 = [t100 * t107 + t101 * t113, -t100 * t105 + t101 * t112, t101 * t109 + t103 * t106, (t103 * pkin(2) + pkin(7) * t114) * t108 + (-pkin(2) * t114 + t103 * pkin(7)) * t106 + t101 * t115 + t103 * pkin(1) + 0; -t103 * t113 + t99 * t107, -t103 * t112 - t99 * t105, t101 * t106 - t103 * t109, (t101 * pkin(2) - pkin(7) * t111) * t108 + (pkin(2) * t111 + t101 * pkin(7)) * t106 - t103 * t115 + t101 * pkin(1) + 0; t104 * t105 + t106 * t112, t104 * t107 - t106 * t113, -t102 * t108, t104 * pkin(6) + qJ(1) + 0 + (pkin(2) * t106 - pkin(7) * t108) * t102; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:43
	% EndTime: 2020-11-04 20:09:43
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (57->32), mult. (136->56), div. (0->0), fcn. (170->10), ass. (0->29)
	t123 = sin(pkin(5));
	t125 = cos(pkin(5));
	t128 = sin(qJ(2));
	t131 = cos(qJ(2));
	t127 = sin(qJ(3));
	t130 = cos(qJ(3));
	t135 = t130 * pkin(3) + t127 * pkin(8) + pkin(2);
	t133 = -pkin(7) * t131 + t135 * t128;
	t134 = t127 * pkin(3) - t130 * pkin(8) + pkin(6);
	t146 = t134 * t123 - t133 * t125;
	t142 = t123 * t130;
	t141 = t123 * t131;
	t140 = t125 * t128;
	t139 = t125 * t131;
	t138 = t127 * t123;
	t137 = t128 * t130;
	t136 = t130 * t131;
	t132 = pkin(7) * t128 + t135 * t131 + pkin(1);
	t129 = cos(qJ(4));
	t126 = sin(qJ(4));
	t124 = cos(pkin(10));
	t122 = sin(pkin(10));
	t121 = t125 * t137 - t138;
	t120 = t123 * t137 + t125 * t127;
	t119 = t122 * t139 + t124 * t128;
	t118 = t122 * t128 - t124 * t139;
	t117 = -t122 * t121 + t124 * t136;
	t116 = t124 * t121 + t122 * t136;
	t1 = [t117 * t129 + t119 * t126, -t117 * t126 + t119 * t129, -(t122 * t140 - t124 * t131) * t127 - t122 * t142, t146 * t122 + t132 * t124 + 0; t116 * t129 + t118 * t126, -t116 * t126 + t118 * t129, (t122 * t131 + t124 * t140) * t127 + t124 * t142, t132 * t122 - t146 * t124 + 0; t120 * t129 - t126 * t141, -t120 * t126 - t129 * t141, -t125 * t130 + t128 * t138, t133 * t123 + t134 * t125 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:43
	% EndTime: 2020-11-04 20:09:43
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (80->45), mult. (150->70), div. (0->0), fcn. (196->12), ass. (0->33)
	t179 = pkin(4) * sin(qJ(4));
	t162 = sin(pkin(5));
	t178 = t162 * pkin(6);
	t161 = sin(pkin(10));
	t164 = cos(pkin(5));
	t177 = t161 * t164;
	t166 = sin(qJ(3));
	t176 = t162 * t166;
	t168 = cos(qJ(3));
	t175 = t162 * t168;
	t169 = cos(qJ(2));
	t174 = t162 * t169;
	t163 = cos(pkin(10));
	t173 = t163 * t164;
	t167 = sin(qJ(2));
	t172 = t164 * t167;
	t171 = t164 * t169;
	t170 = -pkin(9) - pkin(8);
	t160 = qJ(4) + qJ(5);
	t159 = cos(t160);
	t158 = sin(t160);
	t157 = cos(qJ(4)) * pkin(4) + pkin(3);
	t156 = t164 * t166 + t167 * t175;
	t155 = -t164 * t168 + t167 * t176;
	t154 = -t161 * t172 + t163 * t169;
	t153 = t161 * t171 + t163 * t167;
	t152 = t161 * t169 + t163 * t172;
	t151 = t161 * t167 - t163 * t171;
	t150 = -t154 * t166 + t161 * t175;
	t149 = t154 * t168 + t161 * t176;
	t148 = t152 * t168 - t163 * t176;
	t147 = t152 * t166 + t163 * t175;
	t1 = [t149 * t159 + t153 * t158, -t149 * t158 + t153 * t159, -t150, t149 * t157 + t150 * t170 + t153 * t179 + (t163 * pkin(2) + pkin(7) * t177) * t169 + (-pkin(2) * t177 + t163 * pkin(7)) * t167 + t161 * t178 + t163 * pkin(1) + 0; t148 * t159 + t151 * t158, -t148 * t158 + t151 * t159, t147, t148 * t157 - t147 * t170 + t151 * t179 + (t161 * pkin(2) - pkin(7) * t173) * t169 + (pkin(2) * t173 + t161 * pkin(7)) * t167 - t163 * t178 + t161 * pkin(1) + 0; t156 * t159 - t158 * t174, -t156 * t158 - t159 * t174, t155, t164 * pkin(6) - t155 * t170 + t156 * t157 + qJ(1) + 0 + (pkin(2) * t167 + (-pkin(7) - t179) * t169) * t162; 0, 0, 0, 1;];
	Tc_mdh = t1;
end