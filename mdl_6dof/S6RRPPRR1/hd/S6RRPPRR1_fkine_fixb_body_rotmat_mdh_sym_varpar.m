% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:02
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:11
	% EndTime: 2020-11-04 22:02:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:11
	% EndTime: 2020-11-04 22:02:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t69 = cos(qJ(1));
	t68 = sin(qJ(1));
	t1 = [t69, -t68, 0, 0; t68, t69, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:11
	% EndTime: 2020-11-04 22:02:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t73 = cos(qJ(1));
	t72 = cos(qJ(2));
	t71 = sin(qJ(1));
	t70 = sin(qJ(2));
	t1 = [t73 * t72, -t73 * t70, t71, t73 * pkin(1) + t71 * pkin(7) + 0; t71 * t72, -t71 * t70, -t73, t71 * pkin(1) - t73 * pkin(7) + 0; t70, t72, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:11
	% EndTime: 2020-11-04 22:02:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t80 = cos(qJ(1));
	t79 = sin(qJ(1));
	t78 = -qJ(3) - pkin(7);
	t77 = qJ(2) + pkin(10);
	t76 = cos(t77);
	t75 = sin(t77);
	t74 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t80 * t76, -t80 * t75, t79, t80 * t74 - t78 * t79 + 0; t79 * t76, -t79 * t75, -t80, t79 * t74 + t80 * t78 + 0; t75, t76, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:11
	% EndTime: 2020-11-04 22:02:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->17), mult. (23->17), div. (0->0), fcn. (31->8), ass. (0->11)
	t90 = cos(qJ(1));
	t89 = sin(qJ(1));
	t88 = sin(qJ(2));
	t87 = -qJ(3) - pkin(7);
	t86 = cos(pkin(10));
	t85 = sin(pkin(10));
	t84 = qJ(2) + pkin(10);
	t83 = cos(t84);
	t82 = sin(t84);
	t81 = (pkin(3) * t86 + qJ(4) * t85 + pkin(2)) * cos(qJ(2)) + (-t85 * pkin(3) + qJ(4) * t86) * t88 + pkin(1);
	t1 = [t90 * t83, t89, t90 * t82, t81 * t90 - t87 * t89 + 0; t89 * t83, -t90, t89 * t82, t81 * t89 + t90 * t87 + 0; t82, 0, -t83, t88 * pkin(2) + t82 * pkin(3) - t83 * qJ(4) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:11
	% EndTime: 2020-11-04 22:02:11
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (48->21), mult. (35->21), div. (0->0), fcn. (49->10), ass. (0->16)
	t107 = cos(qJ(5));
	t106 = sin(qJ(5));
	t105 = pkin(3) + pkin(4);
	t104 = cos(qJ(1));
	t103 = sin(qJ(1));
	t102 = sin(qJ(2));
	t101 = cos(pkin(10));
	t100 = sin(pkin(10));
	t99 = qJ(2) + pkin(10);
	t98 = qJ(3) + pkin(7) - pkin(8);
	t97 = cos(t99);
	t96 = sin(t99);
	t93 = -t97 * t106 + t96 * t107;
	t92 = -t96 * t106 - t97 * t107;
	t91 = (qJ(4) * t100 + t105 * t101 + pkin(2)) * cos(qJ(2)) + (qJ(4) * t101 - t100 * t105) * t102 + pkin(1);
	t1 = [-t104 * t92, t104 * t93, -t103, t98 * t103 + t91 * t104 + 0; -t103 * t92, t103 * t93, t104, t91 * t103 - t98 * t104 + 0; t93, t92, 0, t102 * pkin(2) - t97 * qJ(4) + t105 * t96 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:11
	% EndTime: 2020-11-04 22:02:12
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (85->34), mult. (79->39), div. (0->0), fcn. (101->14), ass. (0->25)
	t120 = sin(qJ(6));
	t123 = sin(qJ(1));
	t131 = t123 * t120;
	t124 = cos(qJ(6));
	t130 = t123 * t124;
	t126 = cos(qJ(1));
	t129 = t126 * t120;
	t128 = t126 * t124;
	t117 = qJ(2) + pkin(10);
	t127 = pkin(3) + pkin(4);
	t125 = cos(qJ(5));
	t122 = sin(qJ(2));
	t121 = sin(qJ(5));
	t119 = cos(pkin(10));
	t118 = sin(pkin(10));
	t116 = qJ(3) + pkin(7) - pkin(8);
	t115 = -qJ(5) + t117;
	t114 = cos(t117);
	t113 = sin(t117);
	t112 = t119 * pkin(5) - t118 * pkin(9);
	t111 = t118 * pkin(5) + t119 * pkin(9);
	t110 = t113 * t121 + t114 * t125;
	t109 = t113 * t125 - t114 * t121;
	t108 = (qJ(4) * t118 + t111 * t121 + t112 * t125 + t127 * t119 + pkin(2)) * cos(qJ(2)) + (qJ(4) * t119 - t111 * t125 + t112 * t121 - t118 * t127) * t122 + pkin(1);
	t1 = [t110 * t128 - t131, -t110 * t129 - t130, -t126 * t109, t108 * t126 + t116 * t123 + 0; t110 * t130 + t129, -t110 * t131 + t128, -t123 * t109, t108 * t123 - t116 * t126 + 0; t109 * t124, -t109 * t120, t110, pkin(9) * cos(t115) + pkin(5) * sin(t115) + t127 * t113 - t114 * qJ(4) + t122 * pkin(2) + 0 + pkin(6); 0, 0, 0, 1;];
	Tc_mdh = t1;
end