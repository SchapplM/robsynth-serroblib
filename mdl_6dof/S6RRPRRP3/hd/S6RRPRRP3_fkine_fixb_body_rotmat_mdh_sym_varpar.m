% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:13
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:13:49
	% EndTime: 2020-11-04 22:13:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:13:49
	% EndTime: 2020-11-04 22:13:49
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t64 = cos(qJ(1));
	t63 = sin(qJ(1));
	t1 = [t64, -t63, 0, 0; t63, t64, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:13:49
	% EndTime: 2020-11-04 22:13:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t68 = cos(qJ(1));
	t67 = cos(qJ(2));
	t66 = sin(qJ(1));
	t65 = sin(qJ(2));
	t1 = [t68 * t67, -t68 * t65, t66, t68 * pkin(1) + t66 * pkin(7) + 0; t66 * t67, -t66 * t65, -t68, t66 * pkin(1) - t68 * pkin(7) + 0; t65, t67, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:13:49
	% EndTime: 2020-11-04 22:13:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t75 = cos(qJ(1));
	t74 = sin(qJ(1));
	t73 = -qJ(3) - pkin(7);
	t72 = qJ(2) + pkin(10);
	t71 = cos(t72);
	t70 = sin(t72);
	t69 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t75 * t71, -t75 * t70, t74, t75 * t69 - t73 * t74 + 0; t74 * t71, -t74 * t70, -t75, t74 * t69 + t75 * t73 + 0; t70, t71, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:13:49
	% EndTime: 2020-11-04 22:13:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (37->21), mult. (35->25), div. (0->0), fcn. (48->10), ass. (0->17)
	t83 = sin(qJ(4));
	t85 = sin(qJ(1));
	t91 = t85 * t83;
	t86 = cos(qJ(4));
	t90 = t85 * t86;
	t87 = cos(qJ(1));
	t89 = t87 * t83;
	t88 = t87 * t86;
	t84 = sin(qJ(2));
	t82 = -qJ(3) - pkin(7);
	t81 = cos(pkin(10));
	t80 = sin(pkin(10));
	t79 = qJ(2) + pkin(10);
	t78 = cos(t79);
	t77 = sin(t79);
	t76 = (t81 * pkin(3) + t80 * pkin(8) + pkin(2)) * cos(qJ(2)) + (-t80 * pkin(3) + t81 * pkin(8)) * t84 + pkin(1);
	t1 = [t78 * t88 + t91, -t78 * t89 + t90, t87 * t77, t76 * t87 - t82 * t85 + 0; t78 * t90 - t89, -t78 * t91 - t88, t85 * t77, t76 * t85 + t87 * t82 + 0; t77 * t86, -t77 * t83, -t78, t84 * pkin(2) + t77 * pkin(3) - t78 * pkin(8) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:13:49
	% EndTime: 2020-11-04 22:13:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (55->23), mult. (40->24), div. (0->0), fcn. (53->10), ass. (0->17)
	t102 = sin(qJ(1));
	t99 = qJ(4) + qJ(5);
	t96 = sin(t99);
	t110 = t102 * t96;
	t97 = cos(t99);
	t109 = t102 * t97;
	t103 = cos(qJ(1));
	t108 = t103 * t96;
	t107 = t103 * t97;
	t106 = pkin(4) * sin(qJ(4)) + pkin(7) + qJ(3);
	t104 = -pkin(9) - pkin(8);
	t92 = cos(qJ(4)) * pkin(4) + pkin(3);
	t98 = qJ(2) + pkin(10);
	t94 = sin(t98);
	t95 = cos(t98);
	t105 = -t104 * t94 + t92 * t95 + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t95 * t107 + t110, -t95 * t108 + t109, t103 * t94, t106 * t102 + t105 * t103 + 0; t95 * t109 - t108, -t95 * t110 - t107, t102 * t94, t105 * t102 - t106 * t103 + 0; t94 * t97, -t94 * t96, -t95, t94 * t92 + t95 * t104 + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:13:49
	% EndTime: 2020-11-04 22:13:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (68->25), mult. (45->26), div. (0->0), fcn. (58->10), ass. (0->17)
	t120 = qJ(4) + qJ(5);
	t116 = sin(t120);
	t122 = sin(qJ(1));
	t129 = t122 * t116;
	t117 = cos(t120);
	t128 = t122 * t117;
	t123 = cos(qJ(1));
	t127 = t123 * t116;
	t126 = t123 * t117;
	t125 = pkin(5) * t116 + pkin(4) * sin(qJ(4)) + qJ(3) + pkin(7);
	t111 = pkin(5) * t117 + cos(qJ(4)) * pkin(4) + pkin(3);
	t119 = qJ(2) + pkin(10);
	t114 = sin(t119);
	t115 = cos(t119);
	t118 = -qJ(6) - pkin(9) - pkin(8);
	t124 = t111 * t115 - t114 * t118 + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t115 * t126 + t129, -t115 * t127 + t128, t123 * t114, t125 * t122 + t124 * t123 + 0; t115 * t128 - t127, -t115 * t129 - t126, t122 * t114, t124 * t122 - t125 * t123 + 0; t114 * t117, -t114 * t116, -t115, t114 * t111 + t115 * t118 + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end