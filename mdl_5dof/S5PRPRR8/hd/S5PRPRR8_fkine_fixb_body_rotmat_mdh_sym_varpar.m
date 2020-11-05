% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRR8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:01
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:01:24
	% EndTime: 2020-11-04 20:01:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:01:24
	% EndTime: 2020-11-04 20:01:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t80 = cos(pkin(9));
	t79 = sin(pkin(9));
	t1 = [t80, -t79, 0, 0; t79, t80, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:01:24
	% EndTime: 2020-11-04 20:01:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t81 = sin(pkin(9));
	t82 = sin(pkin(5));
	t90 = t81 * t82;
	t83 = cos(pkin(9));
	t89 = t83 * t82;
	t84 = cos(pkin(5));
	t85 = sin(qJ(2));
	t88 = t84 * t85;
	t86 = cos(qJ(2));
	t87 = t84 * t86;
	t1 = [-t81 * t88 + t83 * t86, -t81 * t87 - t83 * t85, t90, t83 * pkin(1) + pkin(6) * t90 + 0; t81 * t86 + t83 * t88, -t81 * t85 + t83 * t87, -t89, t81 * pkin(1) - pkin(6) * t89 + 0; t82 * t85, t82 * t86, t84, t84 * pkin(6) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:01:24
	% EndTime: 2020-11-04 20:01:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (23->23), mult. (43->34), div. (0->0), fcn. (56->6), ass. (0->15)
	t91 = sin(pkin(9));
	t104 = t91 * pkin(2);
	t93 = cos(pkin(9));
	t103 = t93 * pkin(2);
	t92 = sin(pkin(5));
	t102 = t91 * t92;
	t101 = t93 * t92;
	t94 = cos(pkin(5));
	t95 = sin(qJ(2));
	t100 = t94 * t95;
	t96 = cos(qJ(2));
	t99 = t94 * t96;
	t98 = t91 * qJ(3);
	t97 = t93 * qJ(3);
	t1 = [t102, t91 * t100 - t93 * t96, t91 * t99 + t93 * t95, (t94 * t98 + t103) * t96 + (-t94 * t104 + t97) * t95 + pkin(6) * t102 + t93 * pkin(1) + 0; -t101, -t93 * t100 - t91 * t96, t91 * t95 - t93 * t99, (-t94 * t97 + t104) * t96 + (t94 * t103 + t98) * t95 - pkin(6) * t101 + t91 * pkin(1) + 0; t94, -t92 * t95, -t92 * t96, t94 * pkin(6) + qJ(1) + 0 + (pkin(2) * t95 - qJ(3) * t96) * t92; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:01:24
	% EndTime: 2020-11-04 20:01:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (36->28), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->22)
	t107 = sin(pkin(9));
	t125 = t107 * qJ(3);
	t108 = sin(pkin(5));
	t111 = sin(qJ(4));
	t124 = t108 * t111;
	t113 = cos(qJ(4));
	t123 = t108 * t113;
	t114 = cos(qJ(2));
	t122 = t108 * t114;
	t115 = pkin(3) + pkin(6);
	t121 = t108 * t115;
	t109 = cos(pkin(9));
	t120 = t109 * qJ(3);
	t110 = cos(pkin(5));
	t112 = sin(qJ(2));
	t119 = t110 * t112;
	t118 = t110 * t114;
	t116 = pkin(2) + pkin(7);
	t117 = t110 * t116;
	t106 = t107 * t118 + t109 * t112;
	t105 = t107 * t112 - t109 * t118;
	t1 = [t106 * t111 + t107 * t123, t106 * t113 - t107 * t124, -t107 * t119 + t109 * t114, (t109 * t116 + t110 * t125) * t114 + (-t107 * t117 + t120) * t112 + t107 * t121 + t109 * pkin(1) + 0; t105 * t111 - t109 * t123, t105 * t113 + t109 * t124, t107 * t114 + t109 * t119, (t107 * t116 - t110 * t120) * t114 + (t109 * t117 + t125) * t112 - t109 * t121 + t107 * pkin(1) + 0; t110 * t113 - t111 * t122, -t110 * t111 - t113 * t122, t108 * t112, t115 * t110 + qJ(1) + 0 + (-qJ(3) * t114 + t112 * t116) * t108; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:01:24
	% EndTime: 2020-11-04 20:01:24
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (65->34), mult. (140->56), div. (0->0), fcn. (174->10), ass. (0->30)
	t132 = sin(pkin(5));
	t134 = cos(pkin(5));
	t137 = sin(qJ(2));
	t140 = cos(qJ(2));
	t142 = pkin(2) + pkin(7);
	t136 = sin(qJ(4));
	t139 = cos(qJ(4));
	t147 = t136 * pkin(4) - pkin(8) * t139 + qJ(3);
	t144 = -t142 * t137 + t147 * t140;
	t146 = t139 * pkin(4) + t136 * pkin(8) + pkin(3) + pkin(6);
	t158 = t146 * t132 + t144 * t134;
	t155 = t132 * t136;
	t154 = t132 * t137;
	t153 = t134 * t137;
	t152 = t134 * t140;
	t151 = t136 * t137;
	t150 = t136 * t140;
	t149 = t139 * t132;
	t145 = t134 * t150 + t149;
	t143 = t147 * t137 + t142 * t140 + pkin(1);
	t138 = cos(qJ(5));
	t135 = sin(qJ(5));
	t133 = cos(pkin(9));
	t131 = sin(pkin(9));
	t130 = -t132 * t150 + t134 * t139;
	t129 = t131 * t153 - t133 * t140;
	t128 = t131 * t140 + t133 * t153;
	t127 = -t131 * t151 + t145 * t133;
	t126 = t145 * t131 + t133 * t151;
	t1 = [t126 * t138 - t135 * t129, -t126 * t135 - t138 * t129, t131 * t155 - (t131 * t152 + t133 * t137) * t139, t158 * t131 + t143 * t133 + 0; -t127 * t138 + t128 * t135, t127 * t135 + t128 * t138, -t133 * t155 + (-t131 * t137 + t133 * t152) * t139, t143 * t131 - t158 * t133 + 0; t130 * t138 + t135 * t154, -t130 * t135 + t138 * t154, t134 * t136 + t140 * t149, -t144 * t132 + t146 * t134 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end