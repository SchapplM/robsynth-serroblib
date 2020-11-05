% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRPR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:04
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:04:27
	% EndTime: 2020-11-04 20:04:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:04:27
	% EndTime: 2020-11-04 20:04:27
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t83 = cos(pkin(9));
	t82 = sin(pkin(9));
	t1 = [t83, -t82, 0, 0; t82, t83, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:04:27
	% EndTime: 2020-11-04 20:04:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t84 = sin(pkin(9));
	t85 = sin(pkin(5));
	t93 = t84 * t85;
	t86 = cos(pkin(9));
	t92 = t86 * t85;
	t87 = cos(pkin(5));
	t88 = sin(qJ(2));
	t91 = t87 * t88;
	t89 = cos(qJ(2));
	t90 = t87 * t89;
	t1 = [-t84 * t91 + t86 * t89, -t84 * t90 - t86 * t88, t93, t86 * pkin(1) + pkin(6) * t93 + 0; t84 * t89 + t86 * t91, -t84 * t88 + t86 * t90, -t92, t84 * pkin(1) - pkin(6) * t92 + 0; t85 * t88, t85 * t89, t87, t87 * pkin(6) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:04:27
	% EndTime: 2020-11-04 20:04:27
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->46), div. (0->0), fcn. (85->8), ass. (0->20)
	t97 = sin(pkin(5));
	t112 = t97 * pkin(6);
	t96 = sin(pkin(9));
	t99 = cos(pkin(5));
	t111 = t96 * t99;
	t98 = cos(pkin(9));
	t110 = t98 * t99;
	t100 = sin(qJ(3));
	t109 = t100 * t97;
	t102 = cos(qJ(3));
	t108 = t102 * t97;
	t101 = sin(qJ(2));
	t107 = t96 * t101;
	t103 = cos(qJ(2));
	t106 = t96 * t103;
	t105 = t98 * t101;
	t104 = t98 * t103;
	t95 = -t99 * t107 + t104;
	t94 = t99 * t105 + t106;
	t1 = [t95 * t102 + t96 * t109, -t95 * t100 + t96 * t108, t99 * t106 + t105, (t98 * pkin(2) + pkin(7) * t111) * t103 + (-pkin(2) * t111 + t98 * pkin(7)) * t101 + t96 * t112 + t98 * pkin(1) + 0; t94 * t102 - t98 * t109, -t94 * t100 - t98 * t108, -t99 * t104 + t107, (t96 * pkin(2) - pkin(7) * t110) * t103 + (pkin(2) * t110 + t96 * pkin(7)) * t101 - t98 * t112 + t96 * pkin(1) + 0; t99 * t100 + t101 * t108, -t101 * t109 + t99 * t102, -t97 * t103, t99 * pkin(6) + qJ(1) + 0 + (pkin(2) * t101 - pkin(7) * t103) * t97; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:04:27
	% EndTime: 2020-11-04 20:04:27
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (57->39), mult. (134->66), div. (0->0), fcn. (178->10), ass. (0->29)
	t125 = sin(pkin(5));
	t140 = t125 * pkin(6);
	t124 = sin(pkin(9));
	t128 = cos(pkin(5));
	t139 = t124 * t128;
	t129 = sin(qJ(3));
	t138 = t125 * t129;
	t131 = cos(qJ(3));
	t137 = t125 * t131;
	t132 = cos(qJ(2));
	t136 = t125 * t132;
	t127 = cos(pkin(9));
	t135 = t127 * t128;
	t130 = sin(qJ(2));
	t134 = t128 * t130;
	t133 = t128 * t132;
	t126 = cos(pkin(10));
	t123 = sin(pkin(10));
	t122 = t128 * t129 + t130 * t137;
	t121 = -t128 * t131 + t130 * t138;
	t120 = -t124 * t134 + t127 * t132;
	t119 = t124 * t133 + t127 * t130;
	t118 = t124 * t132 + t127 * t134;
	t117 = t124 * t130 - t127 * t133;
	t116 = -t120 * t129 + t124 * t137;
	t115 = t120 * t131 + t124 * t138;
	t114 = t118 * t131 - t127 * t138;
	t113 = t118 * t129 + t127 * t137;
	t1 = [t115 * t126 + t119 * t123, -t115 * t123 + t119 * t126, -t116, t115 * pkin(3) - t116 * qJ(4) + (t127 * pkin(2) + pkin(7) * t139) * t132 + (-pkin(2) * t139 + t127 * pkin(7)) * t130 + t124 * t140 + t127 * pkin(1) + 0; t114 * t126 + t117 * t123, -t114 * t123 + t117 * t126, t113, t114 * pkin(3) + t113 * qJ(4) + (t124 * pkin(2) - pkin(7) * t135) * t132 + (pkin(2) * t135 + t124 * pkin(7)) * t130 - t127 * t140 + t124 * pkin(1) + 0; t122 * t126 - t123 * t136, -t122 * t123 - t126 * t136, t121, t122 * pkin(3) + t128 * pkin(6) + t121 * qJ(4) + qJ(1) + 0 + (pkin(2) * t130 - pkin(7) * t132) * t125; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:04:27
	% EndTime: 2020-11-04 20:04:27
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (90->36), mult. (149->58), div. (0->0), fcn. (183->12), ass. (0->33)
	t153 = sin(pkin(5));
	t155 = cos(pkin(5));
	t147 = sin(pkin(10)) * pkin(4) + pkin(7);
	t158 = sin(qJ(2));
	t160 = cos(qJ(2));
	t148 = cos(pkin(10)) * pkin(4) + pkin(3);
	t156 = qJ(4) + pkin(8);
	t157 = sin(qJ(3));
	t159 = cos(qJ(3));
	t164 = t148 * t159 + t156 * t157 + pkin(2);
	t162 = -t147 * t160 + t164 * t158;
	t163 = t148 * t157 - t156 * t159 + pkin(6);
	t175 = t163 * t153 - t162 * t155;
	t171 = t153 * t159;
	t170 = t153 * t160;
	t169 = t155 * t158;
	t168 = t155 * t160;
	t167 = t157 * t153;
	t166 = t158 * t159;
	t165 = t159 * t160;
	t161 = t147 * t158 + t164 * t160 + pkin(1);
	t154 = cos(pkin(9));
	t152 = sin(pkin(9));
	t151 = pkin(10) + qJ(5);
	t150 = cos(t151);
	t149 = sin(t151);
	t146 = t153 * t166 + t155 * t157;
	t145 = -t155 * t166 + t167;
	t144 = t152 * t168 + t154 * t158;
	t143 = t152 * t158 - t154 * t168;
	t142 = t152 * t145 + t154 * t165;
	t141 = -t154 * t145 + t152 * t165;
	t1 = [t142 * t150 + t144 * t149, -t142 * t149 + t144 * t150, -(t152 * t169 - t154 * t160) * t157 - t152 * t171, t175 * t152 + t161 * t154 + 0; t141 * t150 + t143 * t149, -t141 * t149 + t143 * t150, (t152 * t160 + t154 * t169) * t157 + t154 * t171, t161 * t152 - t175 * t154 + 0; t146 * t150 - t149 * t170, -t146 * t149 - t150 * t170, -t155 * t159 + t158 * t167, t162 * t153 + t163 * t155 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end