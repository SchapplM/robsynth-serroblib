% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRP7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:14
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:58
	% EndTime: 2020-11-04 22:14:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:58
	% EndTime: 2020-11-04 22:14:58
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t65 = cos(qJ(1));
	t64 = sin(qJ(1));
	t1 = [t65, -t64, 0, 0; t64, t65, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:58
	% EndTime: 2020-11-04 22:14:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t69 = cos(qJ(1));
	t68 = cos(qJ(2));
	t67 = sin(qJ(1));
	t66 = sin(qJ(2));
	t1 = [t69 * t68, -t69 * t66, t67, t69 * pkin(1) + t67 * pkin(7) + 0; t67 * t68, -t67 * t66, -t69, t67 * pkin(1) - t69 * pkin(7) + 0; t66, t68, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:58
	% EndTime: 2020-11-04 22:14:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->11), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t74 = cos(qJ(1));
	t73 = cos(qJ(2));
	t72 = sin(qJ(1));
	t71 = sin(qJ(2));
	t70 = t73 * pkin(2) + t71 * qJ(3) + pkin(1);
	t1 = [t74 * t73, t72, t74 * t71, t72 * pkin(7) + t70 * t74 + 0; t72 * t73, -t74, t72 * t71, -t74 * pkin(7) + t70 * t72 + 0; t71, 0, -t73, t71 * pkin(2) - t73 * qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:58
	% EndTime: 2020-11-04 22:14:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (23->16), mult. (26->16), div. (0->0), fcn. (40->6), ass. (0->12)
	t87 = cos(qJ(4));
	t86 = sin(qJ(4));
	t85 = pkin(2) + pkin(3);
	t84 = pkin(7) - pkin(8);
	t83 = cos(qJ(1));
	t82 = cos(qJ(2));
	t81 = sin(qJ(1));
	t80 = sin(qJ(2));
	t77 = t80 * qJ(3) + t85 * t82 + pkin(1);
	t76 = t80 * t87 - t82 * t86;
	t75 = -t80 * t86 - t82 * t87;
	t1 = [-t83 * t75, t83 * t76, -t81, t77 * t83 + t84 * t81 + 0; -t81 * t75, t81 * t76, t83, t77 * t81 - t84 * t83 + 0; t76, t75, 0, -t82 * qJ(3) + t85 * t80 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:58
	% EndTime: 2020-11-04 22:14:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (41->24), mult. (56->28), div. (0->0), fcn. (78->8), ass. (0->19)
	t94 = sin(qJ(5));
	t97 = sin(qJ(1));
	t106 = t97 * t94;
	t98 = cos(qJ(5));
	t105 = t97 * t98;
	t101 = cos(qJ(1));
	t104 = t101 * t94;
	t103 = t101 * t98;
	t102 = pkin(7) - pkin(8);
	t100 = cos(qJ(2));
	t99 = cos(qJ(4));
	t96 = sin(qJ(2));
	t95 = sin(qJ(4));
	t92 = -t95 * pkin(4) + t99 * pkin(9) - qJ(3);
	t91 = t99 * pkin(4) + t95 * pkin(9) + pkin(2) + pkin(3);
	t90 = t100 * t99 + t96 * t95;
	t89 = -t100 * t95 + t96 * t99;
	t88 = t91 * t100 - t92 * t96 + pkin(1);
	t1 = [t90 * t103 - t106, -t90 * t104 - t105, -t101 * t89, t88 * t101 + t102 * t97 + 0; t90 * t105 + t104, -t90 * t106 + t103, -t97 * t89, -t102 * t101 + t88 * t97 + 0; t89 * t98, -t89 * t94, t90, t92 * t100 + t91 * t96 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:14:58
	% EndTime: 2020-11-04 22:14:58
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (57->27), mult. (72->32), div. (0->0), fcn. (94->8), ass. (0->20)
	t114 = sin(qJ(5));
	t117 = sin(qJ(1));
	t126 = t117 * t114;
	t118 = cos(qJ(5));
	t125 = t117 * t118;
	t121 = cos(qJ(1));
	t124 = t121 * t114;
	t123 = t121 * t118;
	t112 = t118 * pkin(5) + qJ(6) * t114 + pkin(4);
	t115 = sin(qJ(4));
	t119 = cos(qJ(4));
	t122 = t119 * pkin(9) - t112 * t115 - qJ(3);
	t120 = cos(qJ(2));
	t116 = sin(qJ(2));
	t111 = -t114 * pkin(5) + qJ(6) * t118 + pkin(7) - pkin(8);
	t110 = t116 * t115 + t120 * t119;
	t109 = -t120 * t115 + t116 * t119;
	t108 = t115 * pkin(9) + t112 * t119 + pkin(2) + pkin(3);
	t107 = t108 * t120 - t122 * t116 + pkin(1);
	t1 = [t110 * t123 - t126, -t121 * t109, t110 * t124 + t125, t107 * t121 + t111 * t117 + 0; t110 * t125 + t124, -t117 * t109, t110 * t126 - t123, t107 * t117 - t111 * t121 + 0; t109 * t118, t110, t109 * t114, t108 * t116 + t122 * t120 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end