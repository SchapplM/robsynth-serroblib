% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:36
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:25
	% EndTime: 2020-11-04 22:36:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:25
	% EndTime: 2020-11-04 22:36:25
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t58 = cos(qJ(1));
	t57 = sin(qJ(1));
	t1 = [t58, -t57, 0, 0; t57, t58, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:25
	% EndTime: 2020-11-04 22:36:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t62 = cos(qJ(1));
	t61 = cos(qJ(2));
	t60 = sin(qJ(1));
	t59 = sin(qJ(2));
	t1 = [t62 * t61, -t62 * t59, t60, t62 * pkin(1) + t60 * pkin(7) + 0; t60 * t61, -t60 * t59, -t62, t60 * pkin(1) - t62 * pkin(7) + 0; t59, t61, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:25
	% EndTime: 2020-11-04 22:36:25
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t66 = sin(qJ(1));
	t68 = cos(qJ(2));
	t72 = t66 * t68;
	t64 = sin(qJ(3));
	t69 = cos(qJ(1));
	t71 = t69 * t64;
	t67 = cos(qJ(3));
	t70 = t69 * t67;
	t65 = sin(qJ(2));
	t63 = t68 * pkin(2) + t65 * pkin(8) + pkin(1);
	t1 = [t66 * t64 + t68 * t70, t66 * t67 - t68 * t71, t69 * t65, t66 * pkin(7) + t63 * t69 + 0; t67 * t72 - t71, -t64 * t72 - t70, t66 * t65, -t69 * pkin(7) + t63 * t66 + 0; t65 * t67, -t65 * t64, -t68, t65 * pkin(2) - t68 * pkin(8) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:25
	% EndTime: 2020-11-04 22:36:25
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (31->23), div. (0->0), fcn. (44->8), ass. (0->15)
	t80 = sin(qJ(1));
	t81 = cos(qJ(2));
	t86 = t80 * t81;
	t78 = qJ(3) + qJ(4);
	t76 = sin(t78);
	t82 = cos(qJ(1));
	t85 = t82 * t76;
	t77 = cos(t78);
	t84 = t82 * t77;
	t83 = pkin(9) + pkin(8);
	t79 = sin(qJ(2));
	t75 = cos(qJ(3)) * pkin(3) + pkin(2);
	t74 = sin(qJ(3)) * pkin(3) + pkin(7);
	t73 = t75 * t81 + t83 * t79 + pkin(1);
	t1 = [t80 * t76 + t81 * t84, t80 * t77 - t81 * t85, t82 * t79, t73 * t82 + t74 * t80 + 0; t77 * t86 - t85, -t76 * t86 - t84, t80 * t79, t73 * t80 - t74 * t82 + 0; t79 * t77, -t79 * t76, -t81, t79 * t75 - t81 * t83 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:25
	% EndTime: 2020-11-04 22:36:25
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (54->24), mult. (51->30), div. (0->0), fcn. (68->8), ass. (0->18)
	t98 = sin(qJ(1));
	t99 = cos(qJ(2));
	t103 = t98 * t99;
	t100 = cos(qJ(1));
	t102 = t100 * t99;
	t101 = pkin(9) + pkin(8);
	t97 = sin(qJ(2));
	t96 = qJ(3) + qJ(4);
	t95 = cos(t96);
	t94 = sin(t96);
	t93 = cos(qJ(3)) * pkin(3) + pkin(2);
	t92 = sin(qJ(3)) * pkin(3) + pkin(7);
	t91 = t101 * t97 + t93 * t99 + pkin(1);
	t90 = t95 * t102 + t98 * t94;
	t89 = t94 * t102 - t98 * t95;
	t88 = -t100 * t94 + t95 * t103;
	t87 = t100 * t95 + t94 * t103;
	t1 = [t90, t100 * t97, t89, t90 * pkin(4) + t89 * qJ(5) + t91 * t100 + t92 * t98 + 0; t88, t98 * t97, t87, t88 * pkin(4) + t87 * qJ(5) - t92 * t100 + t91 * t98 + 0; t97 * t95, -t99, t97 * t94, -t99 * t101 + pkin(6) + 0 + (pkin(4) * t95 + qJ(5) * t94 + t93) * t97; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:25
	% EndTime: 2020-11-04 22:36:25
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (68->25), mult. (57->29), div. (0->0), fcn. (70->10), ass. (0->22)
	t112 = sin(qJ(4));
	t116 = cos(qJ(4));
	t120 = pkin(4) + pkin(5);
	t106 = qJ(5) * t112 + t120 * t116 + pkin(3);
	t107 = qJ(5) * t116 - t112 * t120;
	t113 = sin(qJ(3));
	t117 = cos(qJ(3));
	t127 = -t106 * t113 + t107 * t117 - pkin(7);
	t115 = sin(qJ(1));
	t118 = cos(qJ(2));
	t125 = t115 * t118;
	t111 = qJ(3) + qJ(4);
	t108 = sin(t111);
	t119 = cos(qJ(1));
	t124 = t119 * t108;
	t109 = cos(t111);
	t123 = t119 * t109;
	t121 = t106 * t117 + t107 * t113 + pkin(2);
	t114 = sin(qJ(2));
	t110 = qJ(6) - pkin(8) - pkin(9);
	t104 = -t110 * t114 + t121 * t118 + pkin(1);
	t1 = [t115 * t108 + t118 * t123, -t115 * t109 + t118 * t124, -t119 * t114, t104 * t119 - t127 * t115 + 0; t109 * t125 - t124, t108 * t125 + t123, -t115 * t114, t104 * t115 + t127 * t119 + 0; t114 * t109, t114 * t108, t118, t110 * t118 + t121 * t114 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end