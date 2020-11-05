% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPPR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:23
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:23:23
	% EndTime: 2020-11-04 22:23:23
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:23:23
	% EndTime: 2020-11-04 22:23:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t66 = cos(qJ(1));
	t65 = sin(qJ(1));
	t1 = [t66, -t65, 0, 0; t65, t66, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:23:23
	% EndTime: 2020-11-04 22:23:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t70 = cos(qJ(1));
	t69 = cos(qJ(2));
	t68 = sin(qJ(1));
	t67 = sin(qJ(2));
	t1 = [t70 * t69, -t70 * t67, t68, t70 * pkin(1) + t68 * pkin(7) + 0; t68 * t69, -t68 * t67, -t70, t68 * pkin(1) - t70 * pkin(7) + 0; t67, t69, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:23:23
	% EndTime: 2020-11-04 22:23:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t74 = sin(qJ(1));
	t76 = cos(qJ(2));
	t80 = t74 * t76;
	t72 = sin(qJ(3));
	t77 = cos(qJ(1));
	t79 = t77 * t72;
	t75 = cos(qJ(3));
	t78 = t77 * t75;
	t73 = sin(qJ(2));
	t71 = t76 * pkin(2) + t73 * pkin(8) + pkin(1);
	t1 = [t74 * t72 + t76 * t78, t74 * t75 - t76 * t79, t77 * t73, t74 * pkin(7) + t71 * t77 + 0; t75 * t80 - t79, -t72 * t80 - t78, t74 * t73, -t77 * pkin(7) + t71 * t74 + 0; t73 * t75, -t73 * t72, -t76, t73 * pkin(2) - t76 * pkin(8) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:23:23
	% EndTime: 2020-11-04 22:23:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (31->23), div. (0->0), fcn. (44->8), ass. (0->15)
	t89 = sin(qJ(1));
	t90 = cos(qJ(2));
	t94 = t89 * t90;
	t86 = qJ(3) + pkin(10);
	t84 = sin(t86);
	t91 = cos(qJ(1));
	t93 = t91 * t84;
	t85 = cos(t86);
	t92 = t91 * t85;
	t88 = sin(qJ(2));
	t87 = qJ(4) + pkin(8);
	t83 = cos(qJ(3)) * pkin(3) + pkin(2);
	t82 = sin(qJ(3)) * pkin(3) + pkin(7);
	t81 = t83 * t90 + t87 * t88 + pkin(1);
	t1 = [t89 * t84 + t90 * t92, t89 * t85 - t90 * t93, t91 * t88, t81 * t91 + t82 * t89 + 0; t85 * t94 - t93, -t84 * t94 - t92, t89 * t88, t81 * t89 - t82 * t91 + 0; t88 * t85, -t88 * t84, -t90, t88 * t83 - t90 * t87 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:23:23
	% EndTime: 2020-11-04 22:23:23
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (54->23), mult. (56->29), div. (0->0), fcn. (69->10), ass. (0->21)
	t110 = cos(qJ(1));
	t101 = qJ(3) + pkin(10);
	t99 = sin(t101);
	t114 = t110 * t99;
	t107 = sin(qJ(1));
	t109 = cos(qJ(2));
	t113 = t107 * t109;
	t100 = cos(t101);
	t112 = t110 * t100;
	t105 = sin(qJ(3));
	t108 = cos(qJ(3));
	t102 = sin(pkin(10));
	t103 = cos(pkin(10));
	t97 = pkin(4) * t103 + qJ(5) * t102 + pkin(3);
	t98 = -t102 * pkin(4) + qJ(5) * t103;
	t111 = t97 * t105 - t98 * t108 + pkin(7);
	t106 = sin(qJ(2));
	t104 = qJ(4) + pkin(8);
	t96 = t98 * t105 + t97 * t108 + pkin(2);
	t95 = t104 * t106 + t96 * t109 + pkin(1);
	t1 = [t107 * t99 + t109 * t112, t110 * t106, -t107 * t100 + t109 * t114, t111 * t107 + t95 * t110 + 0; t100 * t113 - t114, t107 * t106, t99 * t113 + t112, t95 * t107 - t111 * t110 + 0; t106 * t100, -t109, t106 * t99, -t109 * t104 + t96 * t106 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:23:23
	% EndTime: 2020-11-04 22:23:23
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (80->31), mult. (80->42), div. (0->0), fcn. (103->12), ass. (0->27)
	t127 = sin(pkin(10));
	t128 = cos(pkin(10));
	t137 = pkin(4) + pkin(5);
	t121 = qJ(5) * t127 + t137 * t128 + pkin(3);
	t122 = qJ(5) * t128 - t127 * t137;
	t130 = sin(qJ(3));
	t134 = cos(qJ(3));
	t142 = -t121 * t130 + t122 * t134 - pkin(7);
	t132 = sin(qJ(1));
	t135 = cos(qJ(2));
	t140 = t132 * t135;
	t136 = cos(qJ(1));
	t139 = t135 * t136;
	t133 = cos(qJ(6));
	t131 = sin(qJ(2));
	t129 = sin(qJ(6));
	t126 = qJ(3) + pkin(10);
	t125 = qJ(4) + pkin(8) - pkin(9);
	t124 = cos(t126);
	t123 = sin(t126);
	t120 = -t132 * t129 + t133 * t139;
	t119 = t129 * t136 + t133 * t140;
	t118 = t129 * t140 - t133 * t136;
	t117 = t129 * t139 + t132 * t133;
	t116 = t121 * t134 + t122 * t130 + pkin(2);
	t115 = t116 * t135 + t125 * t131 + pkin(1);
	t1 = [t123 * t117 + t124 * t120, -t124 * t117 + t120 * t123, -t136 * t131, t115 * t136 - t142 * t132 + 0; t123 * t118 + t119 * t124, -t124 * t118 + t123 * t119, -t132 * t131, t115 * t132 + t142 * t136 + 0; t131 * (t123 * t129 + t124 * t133), t131 * (t123 * t133 - t124 * t129), t135, t116 * t131 - t125 * t135 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end