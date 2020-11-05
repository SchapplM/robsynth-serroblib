% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPP1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:06
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRPP1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:06:15
	% EndTime: 2020-11-04 22:06:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:06:15
	% EndTime: 2020-11-04 22:06:15
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t65 = cos(qJ(1));
	t64 = sin(qJ(1));
	t1 = [t65, -t64, 0, 0; t64, t65, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:06:15
	% EndTime: 2020-11-04 22:06:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t69 = cos(qJ(1));
	t68 = cos(qJ(2));
	t67 = sin(qJ(1));
	t66 = sin(qJ(2));
	t1 = [t69 * t68, -t69 * t66, t67, t69 * pkin(1) + t67 * pkin(7) + 0; t67 * t68, -t67 * t66, -t69, t67 * pkin(1) - t69 * pkin(7) + 0; t66, t68, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:06:15
	% EndTime: 2020-11-04 22:06:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t76 = cos(qJ(1));
	t75 = sin(qJ(1));
	t74 = -qJ(3) - pkin(7);
	t73 = qJ(2) + pkin(9);
	t72 = cos(t73);
	t71 = sin(t73);
	t70 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t76 * t72, -t76 * t71, t75, t76 * t70 - t74 * t75 + 0; t75 * t72, -t75 * t71, -t76, t75 * t70 + t76 * t74 + 0; t71, t72, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:06:15
	% EndTime: 2020-11-04 22:06:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (37->21), mult. (35->25), div. (0->0), fcn. (48->10), ass. (0->17)
	t84 = sin(qJ(4));
	t86 = sin(qJ(1));
	t92 = t86 * t84;
	t87 = cos(qJ(4));
	t91 = t86 * t87;
	t88 = cos(qJ(1));
	t90 = t88 * t84;
	t89 = t88 * t87;
	t85 = sin(qJ(2));
	t83 = -qJ(3) - pkin(7);
	t82 = cos(pkin(9));
	t81 = sin(pkin(9));
	t80 = qJ(2) + pkin(9);
	t79 = cos(t80);
	t78 = sin(t80);
	t77 = (t82 * pkin(3) + t81 * pkin(8) + pkin(2)) * cos(qJ(2)) + (-t81 * pkin(3) + t82 * pkin(8)) * t85 + pkin(1);
	t1 = [t79 * t89 + t92, -t79 * t90 + t91, t88 * t78, t77 * t88 - t83 * t86 + 0; t79 * t91 - t90, -t79 * t92 - t89, t86 * t78, t77 * t86 + t88 * t83 + 0; t78 * t87, -t78 * t84, -t79, t85 * pkin(2) + t78 * pkin(3) - t79 * pkin(8) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:06:15
	% EndTime: 2020-11-04 22:06:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (55->23), mult. (40->26), div. (0->0), fcn. (53->10), ass. (0->15)
	t104 = sin(qJ(1));
	t100 = qJ(2) + pkin(9);
	t98 = cos(t100);
	t109 = t104 * t98;
	t105 = cos(qJ(1));
	t108 = t105 * t98;
	t107 = pkin(4) * sin(qJ(4)) + pkin(7) + qJ(3);
	t101 = -qJ(5) - pkin(8);
	t93 = cos(qJ(4)) * pkin(4) + pkin(3);
	t96 = sin(t100);
	t106 = -t101 * t96 + t93 * t98 + cos(qJ(2)) * pkin(2) + pkin(1);
	t99 = qJ(4) + pkin(10);
	t97 = cos(t99);
	t95 = sin(t99);
	t1 = [t104 * t95 + t97 * t108, t104 * t97 - t95 * t108, t105 * t96, t107 * t104 + t106 * t105 + 0; -t105 * t95 + t97 * t109, -t105 * t97 - t95 * t109, t104 * t96, t106 * t104 - t107 * t105 + 0; t96 * t97, -t96 * t95, -t98, t96 * t93 + t98 * t101 + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:06:15
	% EndTime: 2020-11-04 22:06:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (80->28), mult. (60->30), div. (0->0), fcn. (77->10), ass. (0->21)
	t120 = qJ(4) + pkin(10);
	t116 = sin(t120);
	t125 = sin(qJ(1));
	t132 = t125 * t116;
	t118 = cos(t120);
	t131 = t125 * t118;
	t126 = cos(qJ(1));
	t130 = t126 * t116;
	t129 = t126 * t118;
	t128 = pkin(4) * sin(qJ(4)) + pkin(7) + qJ(3);
	t114 = cos(qJ(4)) * pkin(4) + pkin(3);
	t121 = qJ(2) + pkin(9);
	t117 = sin(t121);
	t119 = cos(t121);
	t122 = -qJ(5) - pkin(8);
	t127 = t114 * t119 - t117 * t122 + cos(qJ(2)) * pkin(2) + pkin(1);
	t113 = t119 * t129 + t132;
	t112 = t119 * t130 - t131;
	t111 = t119 * t131 - t130;
	t110 = t119 * t132 + t129;
	t1 = [t113, t126 * t117, t112, t113 * pkin(5) + t112 * qJ(6) + t128 * t125 + t127 * t126 + 0; t111, t125 * t117, t110, t111 * pkin(5) + t110 * qJ(6) + t127 * t125 - t128 * t126 + 0; t117 * t118, -t119, t117 * t116, sin(qJ(2)) * pkin(2) + t119 * t122 + pkin(6) + 0 + (pkin(5) * t118 + qJ(6) * t116 + t114) * t117; 0, 0, 0, 1;];
	Tc_mdh = t1;
end