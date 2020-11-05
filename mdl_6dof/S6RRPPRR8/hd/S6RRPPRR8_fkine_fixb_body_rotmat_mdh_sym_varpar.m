% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:04
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:04:46
	% EndTime: 2020-11-04 22:04:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:04:46
	% EndTime: 2020-11-04 22:04:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t65 = cos(qJ(1));
	t64 = sin(qJ(1));
	t1 = [t65, -t64, 0, 0; t64, t65, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:04:46
	% EndTime: 2020-11-04 22:04:46
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
	% StartTime: 2020-11-04 22:04:46
	% EndTime: 2020-11-04 22:04:46
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t74 = sin(qJ(1));
	t75 = cos(qJ(2));
	t79 = t74 * t75;
	t71 = sin(pkin(10));
	t76 = cos(qJ(1));
	t78 = t76 * t71;
	t72 = cos(pkin(10));
	t77 = t76 * t72;
	t73 = sin(qJ(2));
	t70 = t75 * pkin(2) + t73 * qJ(3) + pkin(1);
	t1 = [t74 * t71 + t75 * t77, t74 * t72 - t75 * t78, t76 * t73, t74 * pkin(7) + t70 * t76 + 0; t72 * t79 - t78, -t71 * t79 - t77, t74 * t73, -t76 * pkin(7) + t70 * t74 + 0; t73 * t72, -t73 * t71, -t75, t73 * pkin(2) - t75 * qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:04:46
	% EndTime: 2020-11-04 22:04:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (26->18), mult. (36->25), div. (0->0), fcn. (49->6), ass. (0->13)
	t86 = sin(qJ(1));
	t87 = cos(qJ(2));
	t91 = t86 * t87;
	t83 = sin(pkin(10));
	t88 = cos(qJ(1));
	t90 = t88 * t83;
	t84 = cos(pkin(10));
	t89 = t88 * t84;
	t85 = sin(qJ(2));
	t82 = -t83 * pkin(3) + qJ(4) * t84 - pkin(7);
	t81 = t84 * pkin(3) + t83 * qJ(4) + pkin(2);
	t80 = t85 * qJ(3) + t81 * t87 + pkin(1);
	t1 = [t86 * t83 + t87 * t89, t88 * t85, -t86 * t84 + t87 * t90, t80 * t88 - t82 * t86 + 0; t84 * t91 - t90, t86 * t85, t83 * t91 + t89, t80 * t86 + t82 * t88 + 0; t85 * t84, -t87, t85 * t83, -t87 * qJ(3) + t81 * t85 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:04:46
	% EndTime: 2020-11-04 22:04:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (45->23), mult. (56->30), div. (0->0), fcn. (79->8), ass. (0->18)
	t101 = sin(qJ(1));
	t103 = cos(qJ(2));
	t108 = t101 * t103;
	t104 = cos(qJ(1));
	t107 = t103 * t104;
	t105 = pkin(3) + pkin(4);
	t96 = sin(pkin(10));
	t97 = cos(pkin(10));
	t106 = qJ(4) * t97 - t105 * t96 - pkin(7);
	t102 = cos(qJ(5));
	t100 = sin(qJ(2));
	t99 = sin(qJ(5));
	t98 = qJ(3) - pkin(8);
	t95 = t96 * qJ(4) + t105 * t97 + pkin(2);
	t94 = t102 * t97 + t99 * t96;
	t93 = t102 * t96 - t99 * t97;
	t92 = t98 * t100 + t95 * t103 + pkin(1);
	t1 = [t101 * t93 + t94 * t107, -t101 * t94 + t93 * t107, -t104 * t100, -t106 * t101 + t92 * t104 + 0; -t93 * t104 + t94 * t108, t94 * t104 + t93 * t108, -t101 * t100, t92 * t101 + t106 * t104 + 0; t100 * t94, t100 * t93, t103, t95 * t100 - t98 * t103 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:04:46
	% EndTime: 2020-11-04 22:04:46
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (71->33), mult. (80->43), div. (0->0), fcn. (103->10), ass. (0->24)
	t119 = sin(pkin(10));
	t120 = cos(pkin(10));
	t121 = sin(qJ(5));
	t124 = cos(qJ(5));
	t127 = pkin(3) + pkin(4);
	t133 = -qJ(4) * t120 + t127 * t119 + pkin(7) + (t119 * t124 - t120 * t121) * pkin(5);
	t123 = sin(qJ(1));
	t125 = cos(qJ(2));
	t132 = t123 * t125;
	t126 = cos(qJ(1));
	t131 = t126 * t119;
	t130 = t126 * t120;
	t122 = sin(qJ(2));
	t118 = qJ(5) + qJ(6);
	t117 = qJ(3) - pkin(8) - pkin(9);
	t116 = cos(t118);
	t115 = sin(t118);
	t114 = t123 * t119 + t125 * t130;
	t113 = -t123 * t120 + t125 * t131;
	t112 = t120 * t132 - t131;
	t111 = t119 * t132 + t130;
	t110 = t119 * qJ(4) + t127 * t120 + pkin(2) + (t119 * t121 + t120 * t124) * pkin(5);
	t109 = t110 * t125 + t117 * t122 + pkin(1);
	t1 = [t113 * t115 + t114 * t116, t113 * t116 - t114 * t115, -t126 * t122, t109 * t126 + t133 * t123 + 0; t111 * t115 + t112 * t116, t111 * t116 - t112 * t115, -t123 * t122, t109 * t123 - t133 * t126 + 0; t122 * (t115 * t119 + t116 * t120), -t122 * (t115 * t120 - t116 * t119), t125, t110 * t122 - t117 * t125 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end