% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPPR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:58
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPPR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:58:42
	% EndTime: 2020-11-04 21:58:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:58:42
	% EndTime: 2020-11-04 21:58:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t76 = cos(qJ(1));
	t75 = sin(qJ(1));
	t1 = [t76, -t75, 0, 0; t75, t76, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:58:42
	% EndTime: 2020-11-04 21:58:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t80 = cos(qJ(1));
	t79 = cos(qJ(2));
	t78 = sin(qJ(1));
	t77 = sin(qJ(2));
	t1 = [t80 * t79, -t80 * t77, t78, t80 * pkin(1) + t78 * pkin(7) + 0; t78 * t79, -t78 * t77, -t80, t78 * pkin(1) - t80 * pkin(7) + 0; t77, t79, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:58:42
	% EndTime: 2020-11-04 21:58:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t87 = cos(qJ(1));
	t86 = sin(qJ(1));
	t85 = -qJ(3) - pkin(7);
	t84 = qJ(2) + pkin(9);
	t83 = cos(t84);
	t82 = sin(t84);
	t81 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t87 * t83, -t87 * t82, t86, t87 * t81 - t85 * t86 + 0; t86 * t83, -t86 * t82, -t87, t86 * t81 + t87 * t85 + 0; t82, t83, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:58:42
	% EndTime: 2020-11-04 21:58:42
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (37->21), mult. (35->25), div. (0->0), fcn. (48->10), ass. (0->17)
	t92 = sin(pkin(10));
	t98 = sin(qJ(1));
	t103 = t98 * t92;
	t94 = cos(pkin(10));
	t102 = t98 * t94;
	t99 = cos(qJ(1));
	t101 = t99 * t92;
	t100 = t99 * t94;
	t97 = sin(qJ(2));
	t96 = -qJ(3) - pkin(7);
	t95 = cos(pkin(9));
	t93 = sin(pkin(9));
	t91 = qJ(2) + pkin(9);
	t90 = cos(t91);
	t89 = sin(t91);
	t88 = (pkin(3) * t95 + qJ(4) * t93 + pkin(2)) * cos(qJ(2)) + (-t93 * pkin(3) + qJ(4) * t95) * t97 + pkin(1);
	t1 = [t90 * t100 + t103, -t90 * t101 + t102, t99 * t89, t88 * t99 - t96 * t98 + 0; t90 * t102 - t101, -t90 * t103 - t100, t98 * t89, t88 * t98 + t99 * t96 + 0; t89 * t94, -t89 * t92, -t90, t97 * pkin(2) + t89 * pkin(3) - t90 * qJ(4) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:58:42
	% EndTime: 2020-11-04 21:58:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (60->30), mult. (55->35), div. (0->0), fcn. (64->14), ass. (0->20)
	t112 = sin(pkin(10));
	t117 = sin(qJ(1));
	t122 = t117 * t112;
	t114 = cos(pkin(10));
	t121 = t117 * t114;
	t118 = cos(qJ(1));
	t120 = t118 * t112;
	t119 = t118 * t114;
	t111 = qJ(2) + pkin(9);
	t116 = sin(qJ(2));
	t115 = cos(pkin(9));
	t113 = sin(pkin(9));
	t110 = cos(t111);
	t109 = sin(t111);
	t108 = -pkin(10) + t111;
	t107 = pkin(10) + t111;
	t106 = t114 * pkin(4) + qJ(5) * t112 + pkin(3);
	t105 = -t112 * pkin(4) + qJ(5) * t114 - pkin(7) - qJ(3);
	t104 = (qJ(4) * t113 + t106 * t115 + pkin(2)) * cos(qJ(2)) + (qJ(4) * t115 - t113 * t106) * t116 + pkin(1);
	t1 = [t110 * t119 + t122, t118 * t109, t110 * t120 - t121, t104 * t118 - t105 * t117 + 0; t110 * t121 - t120, t117 * t109, t110 * t122 + t119, t104 * t117 + t105 * t118 + 0; t109 * t114, -t110, t109 * t112, t116 * pkin(2) + t109 * pkin(3) - t110 * qJ(4) + pkin(6) + 0 + (-cos(t107) / 0.2e1 + cos(t108) / 0.2e1) * qJ(5) + (sin(t108) / 0.2e1 + sin(t107) / 0.2e1) * pkin(4); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:58:42
	% EndTime: 2020-11-04 21:58:42
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (84->35), mult. (78->39), div. (0->0), fcn. (94->16), ass. (0->26)
	t132 = sin(pkin(10));
	t134 = cos(pkin(10));
	t137 = sin(qJ(6));
	t140 = cos(qJ(6));
	t124 = t140 * t132 - t137 * t134;
	t139 = sin(qJ(1));
	t147 = t139 * t124;
	t125 = t137 * t132 + t140 * t134;
	t146 = t139 * t125;
	t141 = cos(qJ(1));
	t145 = t141 * t124;
	t144 = t141 * t125;
	t131 = qJ(2) + pkin(9);
	t142 = pkin(4) + pkin(5);
	t143 = qJ(5) * t134 - t142 * t132 - pkin(7) - qJ(3);
	t138 = sin(qJ(2));
	t136 = qJ(4) - pkin(8);
	t135 = cos(pkin(9));
	t133 = sin(pkin(9));
	t130 = cos(t131);
	t129 = sin(t131);
	t128 = -pkin(10) + t131;
	t127 = pkin(10) + t131;
	t126 = qJ(5) * t132 + t142 * t134 + pkin(3);
	t123 = (t126 * t135 + t136 * t133 + pkin(2)) * cos(qJ(2)) + (-t133 * t126 + t136 * t135) * t138 + pkin(1);
	t1 = [t130 * t144 + t147, t130 * t145 - t146, -t141 * t129, t123 * t141 - t143 * t139 + 0; t130 * t146 - t145, t130 * t147 + t144, -t139 * t129, t123 * t139 + t143 * t141 + 0; t129 * t125, t129 * t124, t130, -t136 * t130 + t138 * pkin(2) + t129 * pkin(3) + 0 + pkin(6) + (sin(t128) / 0.2e1 + sin(t127) / 0.2e1) * t142 + (cos(t128) / 0.2e1 - cos(t127) / 0.2e1) * qJ(5); 0, 0, 0, 1;];
	Tc_mdh = t1;
end