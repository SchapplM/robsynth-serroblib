% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPP2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:06
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRPP2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPP2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:06:37
	% EndTime: 2020-11-04 22:06:37
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:06:37
	% EndTime: 2020-11-04 22:06:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t61 = cos(qJ(1));
	t60 = sin(qJ(1));
	t1 = [t61, -t60, 0, 0; t60, t61, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:06:37
	% EndTime: 2020-11-04 22:06:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t65 = cos(qJ(1));
	t64 = cos(qJ(2));
	t63 = sin(qJ(1));
	t62 = sin(qJ(2));
	t1 = [t65 * t64, -t65 * t62, t63, t65 * pkin(1) + t63 * pkin(7) + 0; t63 * t64, -t63 * t62, -t65, t63 * pkin(1) - t65 * pkin(7) + 0; t62, t64, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:06:37
	% EndTime: 2020-11-04 22:06:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t72 = cos(qJ(1));
	t71 = sin(qJ(1));
	t70 = -qJ(3) - pkin(7);
	t69 = qJ(2) + pkin(9);
	t68 = cos(t69);
	t67 = sin(t69);
	t66 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t72 * t68, -t72 * t67, t71, t72 * t66 - t70 * t71 + 0; t71 * t68, -t71 * t67, -t72, t71 * t66 + t72 * t70 + 0; t67, t68, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:06:37
	% EndTime: 2020-11-04 22:06:37
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (37->21), mult. (35->25), div. (0->0), fcn. (48->10), ass. (0->17)
	t80 = sin(qJ(4));
	t82 = sin(qJ(1));
	t88 = t82 * t80;
	t83 = cos(qJ(4));
	t87 = t82 * t83;
	t84 = cos(qJ(1));
	t86 = t84 * t80;
	t85 = t84 * t83;
	t81 = sin(qJ(2));
	t79 = -qJ(3) - pkin(7);
	t78 = cos(pkin(9));
	t77 = sin(pkin(9));
	t76 = qJ(2) + pkin(9);
	t75 = cos(t76);
	t74 = sin(t76);
	t73 = (t78 * pkin(3) + t77 * pkin(8) + pkin(2)) * cos(qJ(2)) + (-t77 * pkin(3) + t78 * pkin(8)) * t81 + pkin(1);
	t1 = [t75 * t85 + t88, -t75 * t86 + t87, t84 * t74, t73 * t84 - t79 * t82 + 0; t75 * t87 - t86, -t75 * t88 - t85, t82 * t74, t73 * t82 + t84 * t79 + 0; t74 * t83, -t74 * t80, -t75, t81 * pkin(2) + t74 * pkin(3) - t75 * pkin(8) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:06:37
	% EndTime: 2020-11-04 22:06:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (52->26), mult. (55->31), div. (0->0), fcn. (72->10), ass. (0->21)
	t100 = sin(qJ(4));
	t102 = sin(qJ(1));
	t108 = t102 * t100;
	t103 = cos(qJ(4));
	t107 = t102 * t103;
	t104 = cos(qJ(1));
	t106 = t104 * t100;
	t105 = t104 * t103;
	t101 = sin(qJ(2));
	t99 = -qJ(3) - pkin(7);
	t98 = cos(pkin(9));
	t97 = sin(pkin(9));
	t96 = qJ(2) + pkin(9);
	t95 = cos(t96);
	t94 = sin(t96);
	t93 = t95 * t105 + t108;
	t92 = t95 * t106 - t107;
	t91 = t95 * t107 - t106;
	t90 = t95 * t108 + t105;
	t89 = (t98 * pkin(3) + t97 * pkin(8) + pkin(2)) * cos(qJ(2)) + (-t97 * pkin(3) + t98 * pkin(8)) * t101 + pkin(1);
	t1 = [t93, t104 * t94, t92, t93 * pkin(4) + t92 * qJ(5) - t99 * t102 + t89 * t104 + 0; t91, t102 * t94, t90, t91 * pkin(4) + t90 * qJ(5) + t89 * t102 + t104 * t99 + 0; t94 * t103, -t95, t94 * t100, t101 * pkin(2) - t95 * pkin(8) + pkin(6) + 0 + (pkin(4) * t103 + qJ(5) * t100 + pkin(3)) * t94; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:06:38
	% EndTime: 2020-11-04 22:06:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (74->33), mult. (67->35), div. (0->0), fcn. (73->14), ass. (0->22)
	t118 = sin(qJ(4));
	t120 = sin(qJ(1));
	t129 = t120 * t118;
	t121 = cos(qJ(4));
	t128 = t120 * t121;
	t122 = cos(qJ(1));
	t127 = t122 * t118;
	t126 = t122 * t121;
	t114 = qJ(2) + pkin(9);
	t123 = pkin(4) + pkin(5);
	t125 = qJ(5) * t118 + t121 * t123 + pkin(3);
	t124 = qJ(5) * t121 - t123 * t118 - pkin(7) - qJ(3);
	t119 = sin(qJ(2));
	t117 = qJ(6) - pkin(8);
	t116 = cos(pkin(9));
	t115 = sin(pkin(9));
	t113 = -qJ(4) + t114;
	t112 = qJ(4) + t114;
	t111 = cos(t114);
	t110 = sin(t114);
	t109 = (-t117 * t115 + t125 * t116 + pkin(2)) * cos(qJ(2)) + (-t125 * t115 - t117 * t116) * t119 + pkin(1);
	t1 = [t111 * t126 + t129, t111 * t127 - t128, -t122 * t110, t109 * t122 - t124 * t120 + 0; t111 * t128 - t127, t111 * t129 + t126, -t120 * t110, t109 * t120 + t124 * t122 + 0; t110 * t121, t110 * t118, t111, t117 * t111 + t119 * pkin(2) + t110 * pkin(3) + 0 + pkin(6) + (sin(t113) / 0.2e1 + sin(t112) / 0.2e1) * t123 + (cos(t113) / 0.2e1 - cos(t112) / 0.2e1) * qJ(5); 0, 0, 0, 1;];
	Tc_mdh = t1;
end