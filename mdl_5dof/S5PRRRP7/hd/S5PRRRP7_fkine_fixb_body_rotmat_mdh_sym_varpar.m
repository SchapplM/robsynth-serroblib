% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRP7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:06
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:58
	% EndTime: 2020-11-04 20:06:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:58
	% EndTime: 2020-11-04 20:06:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t87 = cos(pkin(9));
	t86 = sin(pkin(9));
	t1 = [t87, -t86, 0, 0; t86, t87, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:58
	% EndTime: 2020-11-04 20:06:58
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t88 = sin(pkin(9));
	t89 = sin(pkin(5));
	t97 = t88 * t89;
	t90 = cos(pkin(9));
	t96 = t90 * t89;
	t91 = cos(pkin(5));
	t92 = sin(qJ(2));
	t95 = t91 * t92;
	t93 = cos(qJ(2));
	t94 = t91 * t93;
	t1 = [-t88 * t95 + t90 * t93, -t88 * t94 - t90 * t92, t97, t90 * pkin(1) + pkin(6) * t97 + 0; t88 * t93 + t90 * t95, -t88 * t92 + t90 * t94, -t96, t88 * pkin(1) - pkin(6) * t96 + 0; t89 * t92, t89 * t93, t91, t91 * pkin(6) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:58
	% EndTime: 2020-11-04 20:06:58
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t101 = sin(pkin(5));
	t114 = t101 * pkin(6);
	t100 = sin(pkin(9));
	t103 = cos(pkin(5));
	t113 = t100 * t103;
	t104 = sin(qJ(3));
	t112 = t101 * t104;
	t106 = cos(qJ(3));
	t111 = t101 * t106;
	t102 = cos(pkin(9));
	t110 = t102 * t103;
	t105 = sin(qJ(2));
	t109 = t103 * t105;
	t107 = cos(qJ(2));
	t108 = t103 * t107;
	t99 = t100 * t107 + t102 * t109;
	t98 = t100 * t109 - t102 * t107;
	t1 = [t100 * t112 - t98 * t106, t100 * t111 + t98 * t104, t100 * t108 + t102 * t105, (t102 * pkin(2) + pkin(7) * t113) * t107 + (-pkin(2) * t113 + t102 * pkin(7)) * t105 + t100 * t114 + t102 * pkin(1) + 0; -t102 * t112 + t99 * t106, -t102 * t111 - t99 * t104, t100 * t105 - t102 * t108, (t100 * pkin(2) - pkin(7) * t110) * t107 + (pkin(2) * t110 + t100 * pkin(7)) * t105 - t102 * t114 + t100 * pkin(1) + 0; t103 * t104 + t105 * t111, t103 * t106 - t105 * t112, -t101 * t107, t103 * pkin(6) + qJ(1) + 0 + (pkin(2) * t105 - pkin(7) * t107) * t101; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:58
	% EndTime: 2020-11-04 20:06:58
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (57->32), mult. (136->56), div. (0->0), fcn. (170->10), ass. (0->29)
	t122 = sin(pkin(5));
	t124 = cos(pkin(5));
	t127 = sin(qJ(2));
	t130 = cos(qJ(2));
	t126 = sin(qJ(3));
	t129 = cos(qJ(3));
	t134 = t129 * pkin(3) + t126 * pkin(8) + pkin(2);
	t132 = -pkin(7) * t130 + t134 * t127;
	t133 = t126 * pkin(3) - t129 * pkin(8) + pkin(6);
	t145 = t133 * t122 - t132 * t124;
	t141 = t122 * t129;
	t140 = t122 * t130;
	t139 = t124 * t127;
	t138 = t124 * t130;
	t137 = t126 * t122;
	t136 = t127 * t129;
	t135 = t129 * t130;
	t131 = pkin(7) * t127 + t134 * t130 + pkin(1);
	t128 = cos(qJ(4));
	t125 = sin(qJ(4));
	t123 = cos(pkin(9));
	t121 = sin(pkin(9));
	t120 = t124 * t136 - t137;
	t119 = t122 * t136 + t124 * t126;
	t118 = t121 * t138 + t123 * t127;
	t117 = t121 * t127 - t123 * t138;
	t116 = -t121 * t120 + t123 * t135;
	t115 = t123 * t120 + t121 * t135;
	t1 = [t116 * t128 + t118 * t125, -t116 * t125 + t118 * t128, -(t121 * t139 - t123 * t130) * t126 - t121 * t141, t145 * t121 + t131 * t123 + 0; t115 * t128 + t117 * t125, -t115 * t125 + t117 * t128, (t121 * t130 + t123 * t139) * t126 + t123 * t141, t131 * t121 - t145 * t123 + 0; t119 * t128 - t125 * t140, -t119 * t125 - t128 * t140, -t124 * t129 + t127 * t137, t132 * t122 + t133 * t124 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:06:58
	% EndTime: 2020-11-04 20:06:58
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (68->44), mult. (150->70), div. (0->0), fcn. (196->10), ass. (0->33)
	t158 = sin(pkin(5));
	t177 = t158 * pkin(6);
	t157 = sin(pkin(9));
	t159 = cos(pkin(9));
	t164 = sin(qJ(2));
	t160 = cos(pkin(5));
	t167 = cos(qJ(2));
	t168 = t160 * t167;
	t151 = t157 * t164 - t159 * t168;
	t162 = sin(qJ(4));
	t176 = t151 * t162;
	t153 = t157 * t168 + t159 * t164;
	t175 = t153 * t162;
	t174 = t157 * t160;
	t163 = sin(qJ(3));
	t173 = t158 * t163;
	t166 = cos(qJ(3));
	t172 = t158 * t166;
	t171 = t158 * t167;
	t170 = t159 * t160;
	t169 = t160 * t164;
	t165 = cos(qJ(4));
	t161 = -qJ(5) - pkin(8);
	t156 = t165 * pkin(4) + pkin(3);
	t155 = t160 * t163 + t164 * t172;
	t154 = -t160 * t166 + t164 * t173;
	t152 = t157 * t167 + t159 * t169;
	t150 = t157 * t169 - t159 * t167;
	t149 = -t150 * t166 + t157 * t173;
	t148 = t152 * t166 - t159 * t173;
	t147 = t152 * t163 + t159 * t172;
	t146 = t150 * t163 + t157 * t172;
	t1 = [t149 * t165 + t175, -t149 * t162 + t153 * t165, -t146, t149 * t156 + t146 * t161 + pkin(4) * t175 + (t159 * pkin(2) + pkin(7) * t174) * t167 + (-pkin(2) * t174 + t159 * pkin(7)) * t164 + t157 * t177 + t159 * pkin(1) + 0; t148 * t165 + t176, -t148 * t162 + t151 * t165, t147, t148 * t156 - t147 * t161 + pkin(4) * t176 + (t157 * pkin(2) - pkin(7) * t170) * t167 + (pkin(2) * t170 + t157 * pkin(7)) * t164 - t159 * t177 + t157 * pkin(1) + 0; t155 * t165 - t162 * t171, -t155 * t162 - t165 * t171, t154, t160 * pkin(6) - t154 * t161 + t155 * t156 + qJ(1) + 0 + (pkin(2) * t164 + (-pkin(4) * t162 - pkin(7)) * t167) * t158; 0, 0, 0, 1;];
	Tc_mdh = t1;
end