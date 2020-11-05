% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPPR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:24
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPPR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:31
	% EndTime: 2020-11-04 22:24:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:31
	% EndTime: 2020-11-04 22:24:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t61 = cos(qJ(1));
	t60 = sin(qJ(1));
	t1 = [t61, -t60, 0, 0; t60, t61, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:31
	% EndTime: 2020-11-04 22:24:31
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
	% StartTime: 2020-11-04 22:24:31
	% EndTime: 2020-11-04 22:24:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t69 = sin(qJ(1));
	t71 = cos(qJ(2));
	t75 = t69 * t71;
	t67 = sin(qJ(3));
	t72 = cos(qJ(1));
	t74 = t72 * t67;
	t70 = cos(qJ(3));
	t73 = t72 * t70;
	t68 = sin(qJ(2));
	t66 = t71 * pkin(2) + t68 * pkin(8) + pkin(1);
	t1 = [t69 * t67 + t71 * t73, t69 * t70 - t71 * t74, t72 * t68, t69 * pkin(7) + t66 * t72 + 0; t70 * t75 - t74, -t67 * t75 - t73, t69 * t68, -t72 * pkin(7) + t66 * t69 + 0; t68 * t70, -t68 * t67, -t71, t68 * pkin(2) - t71 * pkin(8) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:31
	% EndTime: 2020-11-04 22:24:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (26->18), mult. (36->25), div. (0->0), fcn. (49->6), ass. (0->13)
	t81 = sin(qJ(1));
	t83 = cos(qJ(2));
	t87 = t81 * t83;
	t79 = sin(qJ(3));
	t84 = cos(qJ(1));
	t86 = t84 * t79;
	t82 = cos(qJ(3));
	t85 = t84 * t82;
	t80 = sin(qJ(2));
	t78 = -t79 * pkin(3) + qJ(4) * t82 - pkin(7);
	t77 = t82 * pkin(3) + t79 * qJ(4) + pkin(2);
	t76 = t80 * pkin(8) + t77 * t83 + pkin(1);
	t1 = [t81 * t79 + t83 * t85, t84 * t80, -t81 * t82 + t83 * t86, t76 * t84 - t78 * t81 + 0; t82 * t87 - t86, t81 * t80, t79 * t87 + t85, t76 * t81 + t78 * t84 + 0; t80 * t82, -t83, t80 * t79, -t83 * pkin(8) + t77 * t80 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:31
	% EndTime: 2020-11-04 22:24:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (45->24), mult. (56->29), div. (0->0), fcn. (79->8), ass. (0->19)
	t97 = sin(qJ(1));
	t99 = cos(qJ(2));
	t105 = t97 * t99;
	t100 = cos(qJ(1));
	t92 = sin(pkin(10));
	t93 = cos(pkin(10));
	t95 = sin(qJ(3));
	t98 = cos(qJ(3));
	t89 = t92 * t98 - t93 * t95;
	t104 = t89 * t100;
	t90 = t92 * t95 + t93 * t98;
	t103 = t90 * t100;
	t101 = pkin(3) + pkin(4);
	t102 = qJ(4) * t98 - t101 * t95 - pkin(7);
	t96 = sin(qJ(2));
	t94 = qJ(5) - pkin(8);
	t91 = t95 * qJ(4) + t101 * t98 + pkin(2);
	t88 = t91 * t99 - t94 * t96 + pkin(1);
	t1 = [t99 * t103 - t97 * t89, -t99 * t104 - t97 * t90, -t100 * t96, t88 * t100 - t102 * t97 + 0; t90 * t105 + t104, -t89 * t105 + t103, -t97 * t96, t102 * t100 + t88 * t97 + 0; t96 * t90, -t96 * t89, t99, t91 * t96 + t94 * t99 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:24:31
	% EndTime: 2020-11-04 22:24:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (71->31), mult. (70->39), div. (0->0), fcn. (93->10), ass. (0->23)
	t120 = sin(qJ(1));
	t122 = cos(qJ(2));
	t127 = t120 * t122;
	t118 = sin(qJ(3));
	t123 = cos(qJ(1));
	t126 = t123 * t118;
	t121 = cos(qJ(3));
	t125 = t123 * t121;
	t112 = cos(pkin(10)) * pkin(5) + pkin(3) + pkin(4);
	t113 = sin(pkin(10)) * pkin(5) + qJ(4);
	t124 = t112 * t118 - t113 * t121 + pkin(7);
	t119 = sin(qJ(2));
	t117 = pkin(10) + qJ(6);
	t116 = qJ(5) - pkin(8) + pkin(9);
	t115 = cos(t117);
	t114 = sin(t117);
	t111 = t120 * t118 + t122 * t125;
	t110 = -t120 * t121 + t122 * t126;
	t109 = t121 * t127 - t126;
	t108 = t118 * t127 + t125;
	t107 = t112 * t121 + t113 * t118 + pkin(2);
	t106 = t107 * t122 - t116 * t119 + pkin(1);
	t1 = [t110 * t114 + t111 * t115, t115 * t110 - t111 * t114, -t123 * t119, t106 * t123 + t124 * t120 + 0; t108 * t114 + t109 * t115, t108 * t115 - t109 * t114, -t120 * t119, t106 * t120 - t124 * t123 + 0; t119 * (t114 * t118 + t115 * t121), -t119 * (t114 * t121 - t115 * t118), t122, t107 * t119 + t116 * t122 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end