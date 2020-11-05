% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPP4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:36
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRPP4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:02
	% EndTime: 2020-11-04 22:36:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:02
	% EndTime: 2020-11-04 22:36:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t60 = cos(qJ(1));
	t59 = sin(qJ(1));
	t1 = [t60, -t59, 0, 0; t59, t60, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:02
	% EndTime: 2020-11-04 22:36:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t64 = cos(qJ(1));
	t63 = cos(qJ(2));
	t62 = sin(qJ(1));
	t61 = sin(qJ(2));
	t1 = [t64 * t63, -t64 * t61, t62, t64 * pkin(1) + t62 * pkin(7) + 0; t62 * t63, -t62 * t61, -t64, t62 * pkin(1) - t64 * pkin(7) + 0; t61, t63, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:02
	% EndTime: 2020-11-04 22:36:02
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t68 = sin(qJ(1));
	t70 = cos(qJ(2));
	t74 = t68 * t70;
	t66 = sin(qJ(3));
	t71 = cos(qJ(1));
	t73 = t71 * t66;
	t69 = cos(qJ(3));
	t72 = t71 * t69;
	t67 = sin(qJ(2));
	t65 = t70 * pkin(2) + t67 * pkin(8) + pkin(1);
	t1 = [t68 * t66 + t70 * t72, t68 * t69 - t70 * t73, t71 * t67, t68 * pkin(7) + t65 * t71 + 0; t69 * t74 - t73, -t66 * t74 - t72, t68 * t67, -t71 * pkin(7) + t65 * t68 + 0; t67 * t69, -t67 * t66, -t70, t67 * pkin(2) - t70 * pkin(8) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:02
	% EndTime: 2020-11-04 22:36:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (31->23), div. (0->0), fcn. (44->8), ass. (0->15)
	t82 = sin(qJ(1));
	t83 = cos(qJ(2));
	t88 = t82 * t83;
	t80 = qJ(3) + qJ(4);
	t78 = sin(t80);
	t84 = cos(qJ(1));
	t87 = t84 * t78;
	t79 = cos(t80);
	t86 = t84 * t79;
	t85 = pkin(9) + pkin(8);
	t81 = sin(qJ(2));
	t77 = cos(qJ(3)) * pkin(3) + pkin(2);
	t76 = sin(qJ(3)) * pkin(3) + pkin(7);
	t75 = t77 * t83 + t85 * t81 + pkin(1);
	t1 = [t82 * t78 + t83 * t86, t82 * t79 - t83 * t87, t84 * t81, t75 * t84 + t76 * t82 + 0; t79 * t88 - t87, -t78 * t88 - t86, t82 * t81, t75 * t82 - t76 * t84 + 0; t81 * t79, -t81 * t78, -t83, t81 * t77 - t83 * t85 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:02
	% EndTime: 2020-11-04 22:36:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->22), mult. (42->25), div. (0->0), fcn. (55->10), ass. (0->16)
	t95 = qJ(3) + qJ(4);
	t104 = pkin(7) + pkin(4) * sin(t95) + sin(qJ(3)) * pkin(3);
	t97 = sin(qJ(1));
	t98 = cos(qJ(2));
	t103 = t97 * t98;
	t93 = pkin(10) + t95;
	t91 = sin(t93);
	t99 = cos(qJ(1));
	t102 = t99 * t91;
	t92 = cos(t93);
	t101 = t99 * t92;
	t89 = pkin(4) * cos(t95) + cos(qJ(3)) * pkin(3) + pkin(2);
	t94 = -qJ(5) - pkin(9) - pkin(8);
	t96 = sin(qJ(2));
	t100 = t89 * t98 - t94 * t96 + pkin(1);
	t1 = [t98 * t101 + t97 * t91, -t98 * t102 + t97 * t92, t99 * t96, t100 * t99 + t104 * t97 + 0; t92 * t103 - t102, -t91 * t103 - t101, t97 * t96, t100 * t97 - t104 * t99 + 0; t96 * t92, -t96 * t91, -t98, t96 * t89 + t98 * t94 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:36:02
	% EndTime: 2020-11-04 22:36:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (87->27), mult. (62->31), div. (0->0), fcn. (79->10), ass. (0->20)
	t115 = qJ(3) + qJ(4);
	t124 = pkin(7) + pkin(4) * sin(t115) + sin(qJ(3)) * pkin(3);
	t117 = sin(qJ(1));
	t118 = cos(qJ(2));
	t123 = t117 * t118;
	t113 = pkin(10) + t115;
	t111 = sin(t113);
	t119 = cos(qJ(1));
	t122 = t119 * t111;
	t112 = cos(t113);
	t121 = t119 * t112;
	t109 = pkin(4) * cos(t115) + cos(qJ(3)) * pkin(3) + pkin(2);
	t114 = -qJ(5) - pkin(9) - pkin(8);
	t116 = sin(qJ(2));
	t120 = t109 * t118 - t114 * t116 + pkin(1);
	t108 = t117 * t111 + t118 * t121;
	t107 = -t117 * t112 + t118 * t122;
	t106 = t112 * t123 - t122;
	t105 = t111 * t123 + t121;
	t1 = [t108, t119 * t116, t107, t108 * pkin(5) + t107 * qJ(6) + t124 * t117 + t120 * t119 + 0; t106, t117 * t116, t105, t106 * pkin(5) + t105 * qJ(6) + t120 * t117 - t124 * t119 + 0; t116 * t112, -t118, t116 * t111, t118 * t114 + pkin(6) + 0 + (pkin(5) * t112 + qJ(6) * t111 + t109) * t116; 0, 0, 0, 1;];
	Tc_mdh = t1;
end