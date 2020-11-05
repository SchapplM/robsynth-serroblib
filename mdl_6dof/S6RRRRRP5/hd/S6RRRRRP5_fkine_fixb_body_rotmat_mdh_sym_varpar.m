% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:43
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:43:47
	% EndTime: 2020-11-04 22:43:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:43:47
	% EndTime: 2020-11-04 22:43:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t59 = cos(qJ(1));
	t58 = sin(qJ(1));
	t1 = [t59, -t58, 0, 0; t58, t59, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:43:47
	% EndTime: 2020-11-04 22:43:47
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
	% StartTime: 2020-11-04 22:43:47
	% EndTime: 2020-11-04 22:43:47
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t67 = sin(qJ(1));
	t69 = cos(qJ(2));
	t73 = t67 * t69;
	t65 = sin(qJ(3));
	t70 = cos(qJ(1));
	t72 = t70 * t65;
	t68 = cos(qJ(3));
	t71 = t70 * t68;
	t66 = sin(qJ(2));
	t64 = t69 * pkin(2) + t66 * pkin(8) + pkin(1);
	t1 = [t67 * t65 + t69 * t71, t67 * t68 - t69 * t72, t70 * t66, t67 * pkin(7) + t64 * t70 + 0; t68 * t73 - t72, -t65 * t73 - t71, t67 * t66, -t70 * pkin(7) + t64 * t67 + 0; t66 * t68, -t66 * t65, -t69, t66 * pkin(2) - t69 * pkin(8) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:43:47
	% EndTime: 2020-11-04 22:43:47
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (35->19), mult. (31->23), div. (0->0), fcn. (44->8), ass. (0->15)
	t81 = sin(qJ(1));
	t82 = cos(qJ(2));
	t87 = t81 * t82;
	t79 = qJ(3) + qJ(4);
	t77 = sin(t79);
	t83 = cos(qJ(1));
	t86 = t83 * t77;
	t78 = cos(t79);
	t85 = t83 * t78;
	t84 = pkin(9) + pkin(8);
	t80 = sin(qJ(2));
	t76 = cos(qJ(3)) * pkin(3) + pkin(2);
	t75 = sin(qJ(3)) * pkin(3) + pkin(7);
	t74 = t76 * t82 + t80 * t84 + pkin(1);
	t1 = [t77 * t81 + t82 * t85, t78 * t81 - t82 * t86, t83 * t80, t74 * t83 + t75 * t81 + 0; t78 * t87 - t86, -t77 * t87 - t85, t81 * t80, t74 * t81 - t75 * t83 + 0; t80 * t78, -t80 * t77, -t82, t76 * t80 - t82 * t84 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:43:47
	% EndTime: 2020-11-04 22:43:47
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->22), mult. (42->25), div. (0->0), fcn. (55->10), ass. (0->16)
	t94 = qJ(3) + qJ(4);
	t103 = pkin(7) + pkin(4) * sin(t94) + sin(qJ(3)) * pkin(3);
	t96 = sin(qJ(1));
	t97 = cos(qJ(2));
	t102 = t96 * t97;
	t92 = qJ(5) + t94;
	t90 = sin(t92);
	t98 = cos(qJ(1));
	t101 = t98 * t90;
	t91 = cos(t92);
	t100 = t98 * t91;
	t88 = pkin(4) * cos(t94) + cos(qJ(3)) * pkin(3) + pkin(2);
	t93 = -pkin(10) - pkin(9) - pkin(8);
	t95 = sin(qJ(2));
	t99 = t88 * t97 - t93 * t95 + pkin(1);
	t1 = [t97 * t100 + t96 * t90, -t97 * t101 + t96 * t91, t98 * t95, t103 * t96 + t99 * t98 + 0; t91 * t102 - t101, -t90 * t102 - t100, t96 * t95, -t103 * t98 + t99 * t96 + 0; t95 * t91, -t95 * t90, -t97, t95 * t88 + t97 * t93 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:43:47
	% EndTime: 2020-11-04 22:43:47
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (76->24), mult. (47->27), div. (0->0), fcn. (60->10), ass. (0->16)
	t110 = qJ(3) + qJ(4);
	t109 = qJ(5) + t110;
	t106 = sin(t109);
	t119 = pkin(7) + pkin(5) * t106 + pkin(4) * sin(t110) + sin(qJ(3)) * pkin(3);
	t112 = sin(qJ(1));
	t113 = cos(qJ(2));
	t118 = t112 * t113;
	t114 = cos(qJ(1));
	t117 = t114 * t106;
	t107 = cos(t109);
	t116 = t114 * t107;
	t104 = pkin(5) * t107 + pkin(4) * cos(t110) + cos(qJ(3)) * pkin(3) + pkin(2);
	t108 = -qJ(6) - pkin(10) - pkin(9) - pkin(8);
	t111 = sin(qJ(2));
	t115 = t104 * t113 - t108 * t111 + pkin(1);
	t1 = [t112 * t106 + t113 * t116, t112 * t107 - t113 * t117, t114 * t111, t119 * t112 + t115 * t114 + 0; t107 * t118 - t117, -t106 * t118 - t116, t112 * t111, t115 * t112 - t119 * t114 + 0; t111 * t107, -t111 * t106, -t113, t111 * t104 + t113 * t108 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end