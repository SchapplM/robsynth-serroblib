% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:30
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:30:42
	% EndTime: 2020-11-04 22:30:42
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:30:42
	% EndTime: 2020-11-04 22:30:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t1 = [t57, -t56, 0, 0; t56, t57, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:30:42
	% EndTime: 2020-11-04 22:30:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t61 = cos(qJ(1));
	t60 = cos(qJ(2));
	t59 = sin(qJ(1));
	t58 = sin(qJ(2));
	t1 = [t61 * t60, -t61 * t58, t59, t61 * pkin(1) + t59 * pkin(7) + 0; t59 * t60, -t59 * t58, -t61, t59 * pkin(1) - t61 * pkin(7) + 0; t58, t60, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:30:42
	% EndTime: 2020-11-04 22:30:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t68 = pkin(8) + pkin(7);
	t67 = cos(qJ(1));
	t66 = sin(qJ(1));
	t65 = qJ(2) + qJ(3);
	t64 = cos(t65);
	t63 = sin(t65);
	t62 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t67 * t64, -t67 * t63, t66, t67 * t62 + t68 * t66 + 0; t66 * t64, -t66 * t63, -t67, t66 * t62 - t67 * t68 + 0; t63, t64, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:30:42
	% EndTime: 2020-11-04 22:30:42
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t73 = sin(pkin(11));
	t75 = sin(qJ(1));
	t82 = t75 * t73;
	t74 = cos(pkin(11));
	t81 = t75 * t74;
	t76 = cos(qJ(1));
	t80 = t76 * t73;
	t79 = t76 * t74;
	t72 = qJ(2) + qJ(3);
	t70 = sin(t72);
	t71 = cos(t72);
	t78 = pkin(3) * t71 + qJ(4) * t70 + cos(qJ(2)) * pkin(2) + pkin(1);
	t77 = pkin(8) + pkin(7);
	t1 = [t71 * t79 + t82, -t71 * t80 + t81, t76 * t70, t77 * t75 + t78 * t76 + 0; t71 * t81 - t80, -t71 * t82 - t79, t75 * t70, t78 * t75 - t76 * t77 + 0; t70 * t74, -t70 * t73, -t71, t70 * pkin(3) - t71 * qJ(4) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:30:42
	% EndTime: 2020-11-04 22:30:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (55->23), mult. (40->24), div. (0->0), fcn. (53->10), ass. (0->17)
	t89 = pkin(11) + qJ(5);
	t85 = sin(t89);
	t93 = sin(qJ(1));
	t101 = t93 * t85;
	t86 = cos(t89);
	t100 = t93 * t86;
	t94 = cos(qJ(1));
	t99 = t94 * t85;
	t98 = t94 * t86;
	t97 = sin(pkin(11)) * pkin(4) + pkin(8) + pkin(7);
	t83 = cos(pkin(11)) * pkin(4) + pkin(3);
	t90 = qJ(2) + qJ(3);
	t87 = sin(t90);
	t88 = cos(t90);
	t92 = -pkin(9) - qJ(4);
	t96 = t83 * t88 - t87 * t92 + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t88 * t98 + t101, -t88 * t99 + t100, t94 * t87, t97 * t93 + t96 * t94 + 0; t88 * t100 - t99, -t88 * t101 - t98, t93 * t87, t96 * t93 - t97 * t94 + 0; t87 * t86, -t87 * t85, -t88, t87 * t83 + t88 * t92 + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:30:42
	% EndTime: 2020-11-04 22:30:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (78->26), mult. (45->26), div. (0->0), fcn. (58->12), ass. (0->18)
	t111 = pkin(11) + qJ(5);
	t107 = qJ(6) + t111;
	t104 = sin(t107);
	t113 = sin(qJ(1));
	t121 = t113 * t104;
	t105 = cos(t107);
	t120 = t113 * t105;
	t114 = cos(qJ(1));
	t119 = t114 * t104;
	t118 = t114 * t105;
	t117 = pkin(5) * sin(t111) + sin(pkin(11)) * pkin(4) + pkin(8) + pkin(7);
	t102 = pkin(5) * cos(t111) + cos(pkin(11)) * pkin(4) + pkin(3);
	t112 = qJ(2) + qJ(3);
	t108 = sin(t112);
	t109 = cos(t112);
	t110 = -pkin(10) - pkin(9) - qJ(4);
	t116 = t102 * t109 - t108 * t110 + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t109 * t118 + t121, -t109 * t119 + t120, t114 * t108, t117 * t113 + t116 * t114 + 0; t109 * t120 - t119, -t109 * t121 - t118, t113 * t108, t116 * t113 - t117 * t114 + 0; t108 * t105, -t108 * t104, -t109, t108 * t102 + t109 * t110 + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end