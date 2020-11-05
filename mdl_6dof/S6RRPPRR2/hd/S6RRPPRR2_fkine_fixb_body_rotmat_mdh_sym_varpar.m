% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:02
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:32
	% EndTime: 2020-11-04 22:02:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:32
	% EndTime: 2020-11-04 22:02:32
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t61 = cos(qJ(1));
	t60 = sin(qJ(1));
	t1 = [t61, -t60, 0, 0; t60, t61, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:32
	% EndTime: 2020-11-04 22:02:32
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
	% StartTime: 2020-11-04 22:02:32
	% EndTime: 2020-11-04 22:02:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t72 = cos(qJ(1));
	t71 = sin(qJ(1));
	t70 = -qJ(3) - pkin(7);
	t69 = qJ(2) + pkin(10);
	t68 = cos(t69);
	t67 = sin(t69);
	t66 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t72 * t68, -t72 * t67, t71, t72 * t66 - t70 * t71 + 0; t71 * t68, -t71 * t67, -t72, t71 * t66 + t72 * t70 + 0; t67, t68, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:32
	% EndTime: 2020-11-04 22:02:32
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (37->21), mult. (35->25), div. (0->0), fcn. (48->10), ass. (0->17)
	t77 = sin(pkin(11));
	t83 = sin(qJ(1));
	t88 = t83 * t77;
	t79 = cos(pkin(11));
	t87 = t83 * t79;
	t84 = cos(qJ(1));
	t86 = t84 * t77;
	t85 = t84 * t79;
	t82 = sin(qJ(2));
	t81 = -qJ(3) - pkin(7);
	t80 = cos(pkin(10));
	t78 = sin(pkin(10));
	t76 = qJ(2) + pkin(10);
	t75 = cos(t76);
	t74 = sin(t76);
	t73 = (pkin(3) * t80 + qJ(4) * t78 + pkin(2)) * cos(qJ(2)) + (-t78 * pkin(3) + qJ(4) * t80) * t82 + pkin(1);
	t1 = [t75 * t85 + t88, -t75 * t86 + t87, t84 * t74, t73 * t84 - t81 * t83 + 0; t75 * t87 - t86, -t75 * t88 - t85, t83 * t74, t73 * t83 + t84 * t81 + 0; t74 * t79, -t74 * t77, -t75, t82 * pkin(2) + t74 * pkin(3) - t75 * qJ(4) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:32
	% EndTime: 2020-11-04 22:02:32
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (55->23), mult. (40->26), div. (0->0), fcn. (53->10), ass. (0->15)
	t100 = sin(qJ(1));
	t96 = qJ(2) + pkin(10);
	t94 = cos(t96);
	t105 = t100 * t94;
	t101 = cos(qJ(1));
	t104 = t101 * t94;
	t103 = sin(pkin(11)) * pkin(4) + qJ(3) + pkin(7);
	t89 = cos(pkin(11)) * pkin(4) + pkin(3);
	t92 = sin(t96);
	t99 = -pkin(8) - qJ(4);
	t102 = t89 * t94 - t92 * t99 + cos(qJ(2)) * pkin(2) + pkin(1);
	t95 = pkin(11) + qJ(5);
	t93 = cos(t95);
	t91 = sin(t95);
	t1 = [t100 * t91 + t93 * t104, t100 * t93 - t91 * t104, t101 * t92, t103 * t100 + t102 * t101 + 0; -t101 * t91 + t93 * t105, -t101 * t93 - t91 * t105, t100 * t92, t102 * t100 - t103 * t101 + 0; t92 * t93, -t92 * t91, -t94, t92 * t89 + t94 * t99 + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:32
	% EndTime: 2020-11-04 22:02:32
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (78->26), mult. (45->26), div. (0->0), fcn. (58->12), ass. (0->18)
	t115 = pkin(11) + qJ(5);
	t113 = qJ(6) + t115;
	t108 = sin(t113);
	t118 = sin(qJ(1));
	t125 = t118 * t108;
	t109 = cos(t113);
	t124 = t118 * t109;
	t119 = cos(qJ(1));
	t123 = t119 * t108;
	t122 = t119 * t109;
	t121 = pkin(5) * sin(t115) + sin(pkin(11)) * pkin(4) + qJ(3) + pkin(7);
	t106 = pkin(5) * cos(t115) + cos(pkin(11)) * pkin(4) + pkin(3);
	t116 = qJ(2) + pkin(10);
	t111 = sin(t116);
	t112 = cos(t116);
	t114 = -pkin(9) - pkin(8) - qJ(4);
	t120 = t106 * t112 - t111 * t114 + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t112 * t122 + t125, -t112 * t123 + t124, t119 * t111, t121 * t118 + t120 * t119 + 0; t112 * t124 - t123, -t112 * t125 - t122, t118 * t111, t120 * t118 - t121 * t119 + 0; t111 * t109, -t111 * t108, -t112, t111 * t106 + t112 * t114 + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end