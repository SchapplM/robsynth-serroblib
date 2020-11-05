% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRP3 (for one body)
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

function Tc_mdh = S6RRRRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:43:02
	% EndTime: 2020-11-04 22:43:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:43:02
	% EndTime: 2020-11-04 22:43:03
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t58 = cos(qJ(1));
	t57 = sin(qJ(1));
	t1 = [t58, -t57, 0, 0; t57, t58, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:43:03
	% EndTime: 2020-11-04 22:43:03
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
	% StartTime: 2020-11-04 22:43:03
	% EndTime: 2020-11-04 22:43:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t69 = pkin(8) + pkin(7);
	t68 = cos(qJ(1));
	t67 = sin(qJ(1));
	t66 = qJ(2) + qJ(3);
	t65 = cos(t66);
	t64 = sin(t66);
	t63 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t68 * t65, -t68 * t64, t67, t68 * t63 + t69 * t67 + 0; t67 * t65, -t67 * t64, -t68, t67 * t63 - t68 * t69 + 0; t64, t65, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:43:03
	% EndTime: 2020-11-04 22:43:03
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t74 = sin(qJ(4));
	t75 = sin(qJ(1));
	t83 = t75 * t74;
	t76 = cos(qJ(4));
	t82 = t75 * t76;
	t77 = cos(qJ(1));
	t81 = t77 * t74;
	t80 = t77 * t76;
	t73 = qJ(2) + qJ(3);
	t71 = sin(t73);
	t72 = cos(t73);
	t79 = pkin(3) * t72 + pkin(9) * t71 + cos(qJ(2)) * pkin(2) + pkin(1);
	t78 = pkin(8) + pkin(7);
	t1 = [t72 * t80 + t83, -t72 * t81 + t82, t77 * t71, t78 * t75 + t79 * t77 + 0; t72 * t82 - t81, -t72 * t83 - t80, t75 * t71, t79 * t75 - t77 * t78 + 0; t71 * t76, -t71 * t74, -t72, t71 * pkin(3) - t72 * pkin(9) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:43:03
	% EndTime: 2020-11-04 22:43:03
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (55->23), mult. (40->24), div. (0->0), fcn. (53->10), ass. (0->17)
	t90 = qJ(4) + qJ(5);
	t86 = sin(t90);
	t93 = sin(qJ(1));
	t102 = t93 * t86;
	t88 = cos(t90);
	t101 = t93 * t88;
	t94 = cos(qJ(1));
	t100 = t94 * t86;
	t99 = t94 * t88;
	t98 = pkin(4) * sin(qJ(4)) + pkin(8) + pkin(7);
	t84 = cos(qJ(4)) * pkin(4) + pkin(3);
	t91 = qJ(2) + qJ(3);
	t87 = sin(t91);
	t89 = cos(t91);
	t95 = -pkin(10) - pkin(9);
	t97 = t84 * t89 - t87 * t95 + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t89 * t99 + t102, -t89 * t100 + t101, t94 * t87, t98 * t93 + t97 * t94 + 0; t89 * t101 - t100, -t89 * t102 - t99, t93 * t87, t97 * t93 - t98 * t94 + 0; t87 * t88, -t87 * t86, -t89, t87 * t84 + t89 * t95 + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:43:03
	% EndTime: 2020-11-04 22:43:03
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (68->25), mult. (45->26), div. (0->0), fcn. (58->10), ass. (0->17)
	t111 = qJ(4) + qJ(5);
	t106 = sin(t111);
	t113 = sin(qJ(1));
	t121 = t113 * t106;
	t108 = cos(t111);
	t120 = t113 * t108;
	t114 = cos(qJ(1));
	t119 = t114 * t106;
	t118 = t114 * t108;
	t117 = pkin(5) * t106 + pkin(4) * sin(qJ(4)) + pkin(8) + pkin(7);
	t103 = pkin(5) * t108 + cos(qJ(4)) * pkin(4) + pkin(3);
	t112 = qJ(2) + qJ(3);
	t107 = sin(t112);
	t109 = cos(t112);
	t110 = -qJ(6) - pkin(10) - pkin(9);
	t116 = t103 * t109 - t107 * t110 + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t109 * t118 + t121, -t109 * t119 + t120, t114 * t107, t117 * t113 + t116 * t114 + 0; t109 * t120 - t119, -t109 * t121 - t118, t113 * t107, t116 * t113 - t117 * t114 + 0; t107 * t108, -t107 * t106, -t109, t107 * t103 + t109 * t110 + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end