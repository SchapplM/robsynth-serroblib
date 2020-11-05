% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:10
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:10:11
	% EndTime: 2020-11-04 22:10:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:10:11
	% EndTime: 2020-11-04 22:10:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t67 = cos(qJ(1));
	t66 = sin(qJ(1));
	t1 = [t67, -t66, 0, 0; t66, t67, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:10:11
	% EndTime: 2020-11-04 22:10:11
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t71 = cos(qJ(1));
	t70 = cos(qJ(2));
	t69 = sin(qJ(1));
	t68 = sin(qJ(2));
	t1 = [t71 * t70, -t71 * t68, t69, t71 * pkin(1) + t69 * pkin(7) + 0; t69 * t70, -t69 * t68, -t71, t69 * pkin(1) - t71 * pkin(7) + 0; t68, t70, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:10:11
	% EndTime: 2020-11-04 22:10:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->11), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t76 = cos(qJ(1));
	t75 = cos(qJ(2));
	t74 = sin(qJ(1));
	t73 = sin(qJ(2));
	t72 = t75 * pkin(2) + t73 * qJ(3) + pkin(1);
	t1 = [t76 * t75, t74, t76 * t73, t74 * pkin(7) + t72 * t76 + 0; t74 * t75, -t76, t74 * t73, -t76 * pkin(7) + t72 * t74 + 0; t73, 0, -t75, t73 * pkin(2) - t75 * qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:10:11
	% EndTime: 2020-11-04 22:10:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (23->16), mult. (26->16), div. (0->0), fcn. (40->6), ass. (0->12)
	t89 = cos(qJ(4));
	t88 = sin(qJ(4));
	t87 = pkin(2) + pkin(3);
	t86 = pkin(7) - pkin(8);
	t85 = cos(qJ(1));
	t84 = cos(qJ(2));
	t83 = sin(qJ(1));
	t82 = sin(qJ(2));
	t79 = t82 * qJ(3) + t87 * t84 + pkin(1);
	t78 = t82 * t89 - t84 * t88;
	t77 = -t82 * t88 - t84 * t89;
	t1 = [-t85 * t77, t85 * t78, -t83, t79 * t85 + t86 * t83 + 0; -t83 * t77, t83 * t78, t85, t79 * t83 - t86 * t85 + 0; t78, t77, 0, -t84 * qJ(3) + t87 * t82 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:10:11
	% EndTime: 2020-11-04 22:10:11
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (45->18), mult. (32->18), div. (0->0), fcn. (46->8), ass. (0->14)
	t98 = qJ(4) + pkin(10);
	t103 = sin(t98);
	t102 = cos(qJ(1));
	t101 = cos(qJ(2));
	t100 = sin(qJ(1));
	t99 = sin(qJ(2));
	t97 = qJ(5) - pkin(7) + pkin(8);
	t96 = cos(t98);
	t95 = sin(qJ(4)) * pkin(4) + qJ(3);
	t94 = cos(qJ(4)) * pkin(4) + pkin(2) + pkin(3);
	t92 = -t101 * t103 + t99 * t96;
	t91 = t101 * t96 + t99 * t103;
	t90 = t94 * t101 + t95 * t99 + pkin(1);
	t1 = [t102 * t91, t102 * t92, -t100, -t97 * t100 + t90 * t102 + 0; t100 * t91, t100 * t92, t102, t90 * t100 + t97 * t102 + 0; t92, -t91, 0, -t95 * t101 + t94 * t99 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:10:11
	% EndTime: 2020-11-04 22:10:11
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (77->33), mult. (71->35), div. (0->0), fcn. (93->15), ass. (0->26)
	t129 = pkin(2) + pkin(3);
	t115 = sin(qJ(6));
	t118 = sin(qJ(1));
	t128 = t118 * t115;
	t119 = cos(qJ(6));
	t127 = t118 * t119;
	t122 = cos(qJ(1));
	t126 = t122 * t115;
	t125 = t122 * t119;
	t124 = qJ(4) + pkin(10);
	t123 = sin(t124);
	t121 = cos(qJ(2));
	t120 = cos(qJ(4));
	t117 = sin(qJ(2));
	t116 = sin(qJ(4));
	t114 = cos(pkin(10));
	t113 = sin(pkin(10));
	t112 = qJ(5) - pkin(7) + pkin(8);
	t111 = -qJ(2) + t124;
	t110 = cos(t124);
	t108 = -t113 * pkin(5) + t114 * pkin(9);
	t107 = t114 * pkin(5) + t113 * pkin(9) + pkin(4);
	t106 = t117 * t110 - t121 * t123;
	t105 = t121 * t110 + t117 * t123;
	t104 = (t107 * t120 + t108 * t116 + t129) * t121 + (t107 * t116 - t108 * t120 + qJ(3)) * t117 + pkin(1);
	t1 = [t105 * t125 - t128, -t105 * t126 - t127, -t122 * t106, t104 * t122 - t112 * t118 + 0; t105 * t127 + t126, -t105 * t128 + t125, -t118 * t106, t104 * t118 + t112 * t122 + 0; t106 * t119, -t106 * t115, t105, pkin(9) * cos(t111) - pkin(5) * sin(t111) + pkin(4) * sin(qJ(2) - qJ(4)) + t129 * t117 - t121 * qJ(3) + 0 + pkin(6); 0, 0, 0, 1;];
	Tc_mdh = t1;
end