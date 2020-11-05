% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRR10V2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:50
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRRR10V2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR10V2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:49:59
	% EndTime: 2020-11-04 22:49:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:49:59
	% EndTime: 2020-11-04 22:49:59
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t65 = cos(qJ(1));
	t64 = sin(qJ(1));
	t1 = [t65, -t64, 0, 0; t64, t65, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:49:59
	% EndTime: 2020-11-04 22:49:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->6), mult. (6->6), div. (0->0), fcn. (14->4), ass. (0->5)
	t69 = cos(qJ(1));
	t68 = cos(qJ(2));
	t67 = sin(qJ(1));
	t66 = sin(qJ(2));
	t1 = [t69 * t68, -t69 * t66, t67, t69 * pkin(1) + 0; t67 * t68, -t67 * t66, -t69, t67 * pkin(1) + 0; t66, t68, 0, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:49:59
	% EndTime: 2020-11-04 22:49:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (15->9), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->7)
	t75 = cos(qJ(1));
	t74 = sin(qJ(1));
	t73 = qJ(2) + qJ(3);
	t72 = cos(t73);
	t71 = sin(t73);
	t70 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t75 * t72, -t75 * t71, t74, t75 * t70 + 0; t74 * t72, -t74 * t71, -t75, t74 * t70 + 0; t71, t72, 0, sin(qJ(2)) * pkin(2) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:49:59
	% EndTime: 2020-11-04 22:49:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (33->16), mult. (31->20), div. (0->0), fcn. (44->8), ass. (0->13)
	t80 = sin(qJ(4));
	t81 = sin(qJ(1));
	t88 = t81 * t80;
	t82 = cos(qJ(4));
	t87 = t81 * t82;
	t83 = cos(qJ(1));
	t86 = t83 * t80;
	t85 = t83 * t82;
	t79 = qJ(2) + qJ(3);
	t77 = sin(t79);
	t78 = cos(t79);
	t84 = pkin(3) * t78 + pkin(5) * t77 + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t78 * t85 + t88, -t78 * t86 + t87, t83 * t77, t84 * t83 + 0; t78 * t87 - t86, -t78 * t88 - t85, t81 * t77, t84 * t81 + 0; t77 * t82, -t77 * t80, -t78, t77 * pkin(3) - t78 * pkin(5) + sin(qJ(2)) * pkin(2) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:49:59
	% EndTime: 2020-11-04 22:49:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (45->20), mult. (52->32), div. (0->0), fcn. (73->10), ass. (0->20)
	t94 = qJ(2) + qJ(3);
	t92 = sin(t94);
	t97 = sin(qJ(1));
	t108 = t92 * t97;
	t99 = cos(qJ(4));
	t107 = t92 * t99;
	t96 = sin(qJ(4));
	t106 = t97 * t96;
	t105 = t97 * t99;
	t100 = cos(qJ(1));
	t104 = t100 * t92;
	t103 = t100 * t96;
	t102 = t100 * t99;
	t93 = cos(t94);
	t101 = pkin(3) * t93 + pkin(5) * t92 + cos(qJ(2)) * pkin(2) + pkin(1);
	t98 = cos(qJ(5));
	t95 = sin(qJ(5));
	t90 = t93 * t102 + t106;
	t89 = t93 * t105 - t103;
	t1 = [t95 * t104 + t90 * t98, t98 * t104 - t90 * t95, t93 * t103 - t105, t101 * t100 + 0; t95 * t108 + t89 * t98, t98 * t108 - t89 * t95, t93 * t106 + t102, t101 * t97 + 0; t98 * t107 - t93 * t95, -t95 * t107 - t93 * t98, t92 * t96, t92 * pkin(3) - t93 * pkin(5) + sin(qJ(2)) * pkin(2) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:49:59
	% EndTime: 2020-11-04 22:49:59
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (81->29), mult. (104->47), div. (0->0), fcn. (143->12), ass. (0->31)
	t122 = qJ(2) + qJ(3);
	t120 = sin(t122);
	t125 = sin(qJ(4));
	t139 = t120 * t125;
	t126 = sin(qJ(1));
	t138 = t120 * t126;
	t129 = cos(qJ(4));
	t137 = t120 * t129;
	t130 = cos(qJ(1));
	t136 = t120 * t130;
	t135 = t126 * t125;
	t134 = t126 * t129;
	t133 = t130 * t125;
	t132 = t130 * t129;
	t121 = cos(t122);
	t131 = pkin(3) * t121 + pkin(5) * t120 + cos(qJ(2)) * pkin(2) + pkin(1);
	t128 = cos(qJ(5));
	t127 = cos(qJ(6));
	t124 = sin(qJ(5));
	t123 = sin(qJ(6));
	t118 = t121 * t132 + t135;
	t117 = t121 * t133 - t134;
	t116 = t121 * t134 - t133;
	t115 = t121 * t135 + t132;
	t114 = -t121 * t124 + t128 * t137;
	t113 = t121 * t128 + t124 * t137;
	t112 = t118 * t128 + t124 * t136;
	t111 = t118 * t124 - t128 * t136;
	t110 = t116 * t128 + t124 * t138;
	t109 = t116 * t124 - t128 * t138;
	t1 = [t112 * t127 + t117 * t123, -t112 * t123 + t117 * t127, t111, t111 * pkin(6) + t131 * t130 + 0; t110 * t127 + t115 * t123, -t110 * t123 + t115 * t127, t109, t109 * pkin(6) + t131 * t126 + 0; t114 * t127 + t123 * t139, -t114 * t123 + t127 * t139, t113, t113 * pkin(6) + t120 * pkin(3) - t121 * pkin(5) + sin(qJ(2)) * pkin(2) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end