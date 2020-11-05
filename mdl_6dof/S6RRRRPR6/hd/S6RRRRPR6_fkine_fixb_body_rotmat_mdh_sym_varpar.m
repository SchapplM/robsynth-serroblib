% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:39
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:39:03
	% EndTime: 2020-11-04 22:39:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:39:03
	% EndTime: 2020-11-04 22:39:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t58 = cos(qJ(1));
	t57 = sin(qJ(1));
	t1 = [t58, -t57, 0, 0; t57, t58, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:39:03
	% EndTime: 2020-11-04 22:39:03
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
	% StartTime: 2020-11-04 22:39:03
	% EndTime: 2020-11-04 22:39:03
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
	% StartTime: 2020-11-04 22:39:03
	% EndTime: 2020-11-04 22:39:03
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
	% StartTime: 2020-11-04 22:39:03
	% EndTime: 2020-11-04 22:39:03
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->22), mult. (42->25), div. (0->0), fcn. (55->10), ass. (0->16)
	t93 = qJ(3) + qJ(4);
	t102 = pkin(7) + pkin(4) * sin(t93) + sin(qJ(3)) * pkin(3);
	t95 = sin(qJ(1));
	t96 = cos(qJ(2));
	t101 = t95 * t96;
	t91 = pkin(11) + t93;
	t89 = sin(t91);
	t97 = cos(qJ(1));
	t100 = t97 * t89;
	t90 = cos(t91);
	t99 = t97 * t90;
	t87 = pkin(4) * cos(t93) + cos(qJ(3)) * pkin(3) + pkin(2);
	t92 = -qJ(5) - pkin(9) - pkin(8);
	t94 = sin(qJ(2));
	t98 = t87 * t96 - t92 * t94 + pkin(1);
	t1 = [t95 * t89 + t96 * t99, -t96 * t100 + t95 * t90, t97 * t94, t102 * t95 + t98 * t97 + 0; t90 * t101 - t100, -t89 * t101 - t99, t95 * t94, -t102 * t97 + t98 * t95 + 0; t94 * t90, -t94 * t89, -t96, t94 * t87 + t96 * t92 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:39:03
	% EndTime: 2020-11-04 22:39:03
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (86->25), mult. (47->27), div. (0->0), fcn. (60->12), ass. (0->17)
	t110 = qJ(3) + qJ(4);
	t108 = pkin(11) + t110;
	t119 = pkin(7) + pkin(5) * sin(t108) + pkin(4) * sin(t110) + sin(qJ(3)) * pkin(3);
	t112 = sin(qJ(1));
	t113 = cos(qJ(2));
	t118 = t112 * t113;
	t107 = qJ(6) + t108;
	t105 = sin(t107);
	t114 = cos(qJ(1));
	t117 = t114 * t105;
	t106 = cos(t107);
	t116 = t114 * t106;
	t103 = pkin(5) * cos(t108) + pkin(4) * cos(t110) + cos(qJ(3)) * pkin(3) + pkin(2);
	t109 = -pkin(10) - qJ(5) - pkin(9) - pkin(8);
	t111 = sin(qJ(2));
	t115 = t103 * t113 - t109 * t111 + pkin(1);
	t1 = [t112 * t105 + t113 * t116, t112 * t106 - t113 * t117, t114 * t111, t119 * t112 + t115 * t114 + 0; t106 * t118 - t117, -t105 * t118 - t116, t112 * t111, t115 * t112 - t119 * t114 + 0; t111 * t106, -t111 * t105, -t113, t111 * t103 + t113 * t109 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end