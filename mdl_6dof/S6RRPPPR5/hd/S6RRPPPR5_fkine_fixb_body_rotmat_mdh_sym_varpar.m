% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPPR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:00
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:00:04
	% EndTime: 2020-11-04 22:00:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:00:04
	% EndTime: 2020-11-04 22:00:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t59 = cos(qJ(1));
	t58 = sin(qJ(1));
	t1 = [t59, -t58, 0, 0; t58, t59, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:00:04
	% EndTime: 2020-11-04 22:00:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t63 = cos(qJ(1));
	t62 = cos(qJ(2));
	t61 = sin(qJ(1));
	t60 = sin(qJ(2));
	t1 = [t63 * t62, -t63 * t60, t61, t63 * pkin(1) + t61 * pkin(7) + 0; t61 * t62, -t61 * t60, -t63, t61 * pkin(1) - t63 * pkin(7) + 0; t60, t62, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:00:04
	% EndTime: 2020-11-04 22:00:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t68 = sin(qJ(1));
	t69 = cos(qJ(2));
	t73 = t68 * t69;
	t65 = sin(pkin(9));
	t70 = cos(qJ(1));
	t72 = t70 * t65;
	t66 = cos(pkin(9));
	t71 = t70 * t66;
	t67 = sin(qJ(2));
	t64 = t69 * pkin(2) + t67 * qJ(3) + pkin(1);
	t1 = [t68 * t65 + t69 * t71, t68 * t66 - t69 * t72, t70 * t67, t68 * pkin(7) + t64 * t70 + 0; t66 * t73 - t72, -t65 * t73 - t71, t68 * t67, -t70 * pkin(7) + t64 * t68 + 0; t67 * t66, -t67 * t65, -t69, t67 * pkin(2) - t69 * qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:00:04
	% EndTime: 2020-11-04 22:00:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (27->19), mult. (36->25), div. (0->0), fcn. (49->6), ass. (0->13)
	t80 = sin(qJ(1));
	t81 = cos(qJ(2));
	t85 = t80 * t81;
	t77 = sin(pkin(9));
	t82 = cos(qJ(1));
	t84 = t82 * t77;
	t78 = cos(pkin(9));
	t83 = t82 * t78;
	t79 = sin(qJ(2));
	t76 = -t77 * pkin(3) + qJ(4) * t78 - pkin(7);
	t75 = t78 * pkin(3) + t77 * qJ(4) + pkin(2);
	t74 = t79 * qJ(3) + t75 * t81 + pkin(1);
	t1 = [t82 * t79, -t80 * t77 - t81 * t83, -t80 * t78 + t81 * t84, t74 * t82 - t76 * t80 + 0; t80 * t79, -t78 * t85 + t84, t77 * t85 + t83, t74 * t80 + t76 * t82 + 0; -t81, -t79 * t78, t79 * t77, -t81 * qJ(3) + t75 * t79 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:00:04
	% EndTime: 2020-11-04 22:00:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->21), mult. (36->25), div. (0->0), fcn. (49->6), ass. (0->15)
	t93 = sin(qJ(1));
	t94 = cos(qJ(2));
	t99 = t93 * t94;
	t88 = sin(pkin(9));
	t95 = cos(qJ(1));
	t98 = t95 * t88;
	t89 = cos(pkin(9));
	t97 = t95 * t89;
	t90 = qJ(5) + pkin(3);
	t96 = qJ(4) * t89 - t90 * t88 - pkin(7);
	t92 = sin(qJ(2));
	t91 = qJ(3) + pkin(4);
	t87 = t88 * qJ(4) + t90 * t89 + pkin(2);
	t86 = t87 * t94 + t91 * t92 + pkin(1);
	t1 = [-t93 * t89 + t94 * t98, -t95 * t92, t93 * t88 + t94 * t97, t86 * t95 - t96 * t93 + 0; t88 * t99 + t97, -t93 * t92, t89 * t99 - t98, t86 * t93 + t96 * t95 + 0; t92 * t88, t94, t92 * t89, t87 * t92 - t91 * t94 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:00:04
	% EndTime: 2020-11-04 22:00:04
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (52->23), mult. (56->30), div. (0->0), fcn. (79->8), ass. (0->19)
	t111 = sin(qJ(1));
	t113 = cos(qJ(2));
	t117 = t111 * t113;
	t114 = cos(qJ(1));
	t116 = t113 * t114;
	t105 = sin(pkin(9));
	t106 = cos(pkin(9));
	t107 = qJ(5) + pkin(3);
	t108 = qJ(4) + pkin(5);
	t115 = t107 * t105 - t108 * t106 + pkin(7);
	t112 = cos(qJ(6));
	t110 = sin(qJ(2));
	t109 = sin(qJ(6));
	t104 = qJ(3) + pkin(4) + pkin(8);
	t103 = -t109 * t105 + t112 * t106;
	t102 = t112 * t105 + t109 * t106;
	t101 = t108 * t105 + t107 * t106 + pkin(2);
	t100 = t101 * t113 + t104 * t110 + pkin(1);
	t1 = [t102 * t116 - t111 * t103, t111 * t102 + t103 * t116, t114 * t110, t100 * t114 + t115 * t111 + 0; t102 * t117 + t103 * t114, -t102 * t114 + t103 * t117, t111 * t110, t100 * t111 - t115 * t114 + 0; t110 * t102, t110 * t103, -t113, t101 * t110 - t104 * t113 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end