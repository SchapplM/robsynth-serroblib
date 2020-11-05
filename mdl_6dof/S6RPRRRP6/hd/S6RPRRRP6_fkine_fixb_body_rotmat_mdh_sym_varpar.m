% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRP6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:52
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:52:41
	% EndTime: 2020-11-04 21:52:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:52:41
	% EndTime: 2020-11-04 21:52:41
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t60 = cos(qJ(1));
	t59 = sin(qJ(1));
	t1 = [t60, -t59, 0, 0; t59, t60, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:52:41
	% EndTime: 2020-11-04 21:52:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t64 = cos(qJ(1));
	t63 = sin(qJ(1));
	t62 = cos(pkin(10));
	t61 = sin(pkin(10));
	t1 = [t64 * t62, -t64 * t61, t63, t64 * pkin(1) + t63 * qJ(2) + 0; t63 * t62, -t63 * t61, -t64, t63 * pkin(1) - t64 * qJ(2) + 0; t61, t62, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:52:41
	% EndTime: 2020-11-04 21:52:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t71 = cos(qJ(1));
	t70 = sin(qJ(1));
	t69 = pkin(7) + qJ(2);
	t68 = pkin(10) + qJ(3);
	t67 = cos(t68);
	t66 = sin(t68);
	t65 = cos(pkin(10)) * pkin(2) + pkin(1);
	t1 = [t71 * t67, -t71 * t66, t70, t71 * t65 + t69 * t70 + 0; t70 * t67, -t70 * t66, -t71, t70 * t65 - t71 * t69 + 0; t66, t67, 0, sin(pkin(10)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:52:41
	% EndTime: 2020-11-04 21:52:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t77 = sin(qJ(4));
	t78 = sin(qJ(1));
	t85 = t78 * t77;
	t79 = cos(qJ(4));
	t84 = t78 * t79;
	t80 = cos(qJ(1));
	t83 = t80 * t77;
	t82 = t80 * t79;
	t75 = pkin(10) + qJ(3);
	t73 = sin(t75);
	t74 = cos(t75);
	t81 = pkin(3) * t74 + pkin(8) * t73 + cos(pkin(10)) * pkin(2) + pkin(1);
	t76 = pkin(7) + qJ(2);
	t1 = [t74 * t82 + t85, -t74 * t83 + t84, t80 * t73, t76 * t78 + t81 * t80 + 0; t74 * t84 - t83, -t74 * t85 - t82, t78 * t73, -t80 * t76 + t81 * t78 + 0; t73 * t79, -t73 * t77, -t74, t73 * pkin(3) - t74 * pkin(8) + sin(pkin(10)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:52:41
	% EndTime: 2020-11-04 21:52:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (55->23), mult. (40->24), div. (0->0), fcn. (53->10), ass. (0->17)
	t93 = qJ(4) + qJ(5);
	t90 = sin(t93);
	t96 = sin(qJ(1));
	t104 = t96 * t90;
	t91 = cos(t93);
	t103 = t96 * t91;
	t97 = cos(qJ(1));
	t102 = t97 * t90;
	t101 = t97 * t91;
	t100 = pkin(4) * sin(qJ(4)) + pkin(7) + qJ(2);
	t87 = cos(qJ(4)) * pkin(4) + pkin(3);
	t92 = pkin(10) + qJ(3);
	t88 = sin(t92);
	t89 = cos(t92);
	t98 = -pkin(9) - pkin(8);
	t99 = t87 * t89 - t88 * t98 + cos(pkin(10)) * pkin(2) + pkin(1);
	t1 = [t89 * t101 + t104, -t89 * t102 + t103, t97 * t88, t100 * t96 + t99 * t97 + 0; t89 * t103 - t102, -t89 * t104 - t101, t96 * t88, -t100 * t97 + t99 * t96 + 0; t88 * t91, -t88 * t90, -t89, t88 * t87 + t89 * t98 + sin(pkin(10)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:52:41
	% EndTime: 2020-11-04 21:52:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (68->25), mult. (45->26), div. (0->0), fcn. (58->10), ass. (0->17)
	t114 = qJ(4) + qJ(5);
	t110 = sin(t114);
	t116 = sin(qJ(1));
	t123 = t116 * t110;
	t111 = cos(t114);
	t122 = t116 * t111;
	t117 = cos(qJ(1));
	t121 = t117 * t110;
	t120 = t117 * t111;
	t119 = pkin(5) * t110 + pkin(4) * sin(qJ(4)) + pkin(7) + qJ(2);
	t105 = pkin(5) * t111 + cos(qJ(4)) * pkin(4) + pkin(3);
	t113 = pkin(10) + qJ(3);
	t108 = sin(t113);
	t109 = cos(t113);
	t112 = -qJ(6) - pkin(9) - pkin(8);
	t118 = t105 * t109 - t108 * t112 + cos(pkin(10)) * pkin(2) + pkin(1);
	t1 = [t109 * t120 + t123, -t109 * t121 + t122, t117 * t108, t119 * t116 + t118 * t117 + 0; t109 * t122 - t121, -t109 * t123 - t120, t116 * t108, t118 * t116 - t119 * t117 + 0; t108 * t111, -t108 * t110, -t109, t108 * t105 + t109 * t112 + sin(pkin(10)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end