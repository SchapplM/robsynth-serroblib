% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRP8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:28
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRP8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:28:32
	% EndTime: 2020-11-04 22:28:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:28:32
	% EndTime: 2020-11-04 22:28:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t61 = cos(qJ(1));
	t60 = sin(qJ(1));
	t1 = [t61, -t60, 0, 0; t60, t61, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:28:32
	% EndTime: 2020-11-04 22:28:32
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
	% StartTime: 2020-11-04 22:28:32
	% EndTime: 2020-11-04 22:28:32
	% DurationCPUTime: 0.03s
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
	% StartTime: 2020-11-04 22:28:32
	% EndTime: 2020-11-04 22:28:32
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
	% StartTime: 2020-11-04 22:28:32
	% EndTime: 2020-11-04 22:28:32
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (45->23), mult. (56->30), div. (0->0), fcn. (79->8), ass. (0->18)
	t95 = sin(qJ(1));
	t98 = cos(qJ(2));
	t104 = t95 * t98;
	t99 = cos(qJ(1));
	t103 = t98 * t99;
	t101 = pkin(3) + pkin(4);
	t93 = sin(qJ(3));
	t97 = cos(qJ(3));
	t102 = qJ(4) * t97 - t101 * t93 - pkin(7);
	t100 = pkin(8) - pkin(9);
	t96 = cos(qJ(5));
	t94 = sin(qJ(2));
	t92 = sin(qJ(5));
	t91 = t93 * qJ(4) + t101 * t97 + pkin(2);
	t90 = t92 * t93 + t96 * t97;
	t89 = -t92 * t97 + t96 * t93;
	t88 = t100 * t94 + t91 * t98 + pkin(1);
	t1 = [t90 * t103 + t95 * t89, t89 * t103 - t95 * t90, -t99 * t94, -t102 * t95 + t88 * t99 + 0; t90 * t104 - t89 * t99, t89 * t104 + t90 * t99, -t95 * t94, t102 * t99 + t88 * t95 + 0; t94 * t90, t94 * t89, t98, -t100 * t98 + t91 * t94 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:28:32
	% EndTime: 2020-11-04 22:28:32
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->25), mult. (66->32), div. (0->0), fcn. (89->8), ass. (0->19)
	t115 = sin(qJ(1));
	t118 = cos(qJ(2));
	t122 = t115 * t118;
	t119 = cos(qJ(1));
	t121 = t118 * t119;
	t116 = cos(qJ(5));
	t109 = t116 * pkin(5) + pkin(3) + pkin(4);
	t112 = sin(qJ(5));
	t110 = t112 * pkin(5) + qJ(4);
	t113 = sin(qJ(3));
	t117 = cos(qJ(3));
	t120 = t109 * t113 - t110 * t117 + pkin(7);
	t114 = sin(qJ(2));
	t111 = qJ(6) - pkin(8) + pkin(9);
	t108 = t112 * t113 + t116 * t117;
	t107 = -t112 * t117 + t116 * t113;
	t106 = t109 * t117 + t110 * t113 + pkin(2);
	t105 = t106 * t118 - t111 * t114 + pkin(1);
	t1 = [t115 * t107 + t108 * t121, t107 * t121 - t115 * t108, -t119 * t114, t105 * t119 + t120 * t115 + 0; -t107 * t119 + t108 * t122, t107 * t122 + t108 * t119, -t115 * t114, t105 * t115 - t120 * t119 + 0; t114 * t108, t114 * t107, t118, t106 * t114 + t111 * t118 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end