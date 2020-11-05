% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:27
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:27
	% EndTime: 2020-11-04 22:27:27
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:27
	% EndTime: 2020-11-04 22:27:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t59 = cos(qJ(1));
	t58 = sin(qJ(1));
	t1 = [t59, -t58, 0, 0; t58, t59, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:27
	% EndTime: 2020-11-04 22:27:27
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
	% StartTime: 2020-11-04 22:27:27
	% EndTime: 2020-11-04 22:27:27
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
	% StartTime: 2020-11-04 22:27:27
	% EndTime: 2020-11-04 22:27:27
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (31->23), div. (0->0), fcn. (44->8), ass. (0->15)
	t82 = sin(qJ(1));
	t83 = cos(qJ(2));
	t87 = t82 * t83;
	t79 = qJ(3) + pkin(10);
	t77 = sin(t79);
	t84 = cos(qJ(1));
	t86 = t84 * t77;
	t78 = cos(t79);
	t85 = t84 * t78;
	t81 = sin(qJ(2));
	t80 = qJ(4) + pkin(8);
	t76 = cos(qJ(3)) * pkin(3) + pkin(2);
	t75 = sin(qJ(3)) * pkin(3) + pkin(7);
	t74 = t76 * t83 + t80 * t81 + pkin(1);
	t1 = [t82 * t77 + t83 * t85, t82 * t78 - t83 * t86, t84 * t81, t74 * t84 + t75 * t82 + 0; t78 * t87 - t86, -t77 * t87 - t85, t82 * t81, t74 * t82 - t75 * t84 + 0; t81 * t78, -t81 * t77, -t83, t81 * t76 - t83 * t80 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:27
	% EndTime: 2020-11-04 22:27:27
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->25), mult. (44->30), div. (0->0), fcn. (57->11), ass. (0->18)
	t105 = sin(pkin(10)) * pkin(4);
	t97 = sin(qJ(1));
	t99 = cos(qJ(2));
	t104 = t97 * t99;
	t100 = cos(qJ(1));
	t103 = t100 * t99;
	t102 = qJ(3) + pkin(10);
	t91 = cos(pkin(10)) * pkin(4) + pkin(3);
	t95 = sin(qJ(3));
	t98 = cos(qJ(3));
	t101 = t98 * t105 + t95 * t91 + pkin(7);
	t96 = sin(qJ(2));
	t93 = pkin(9) + qJ(4) + pkin(8);
	t92 = qJ(5) + t102;
	t90 = cos(t92);
	t89 = sin(t92);
	t88 = (-t95 * t105 + t98 * t91 + pkin(2)) * t99 + t93 * t96 + pkin(1);
	t1 = [t90 * t103 + t97 * t89, -t89 * t103 + t97 * t90, t100 * t96, t88 * t100 + t101 * t97 + 0; -t100 * t89 + t90 * t104, -t100 * t90 - t89 * t104, t97 * t96, -t101 * t100 + t88 * t97 + 0; t96 * t90, -t96 * t89, -t99, t96 * (pkin(4) * cos(t102) + t98 * pkin(3) + pkin(2)) - t99 * t93 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:27:27
	% EndTime: 2020-11-04 22:27:27
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (87->30), mult. (64->35), div. (0->0), fcn. (81->11), ass. (0->23)
	t128 = sin(pkin(10)) * pkin(4);
	t119 = sin(qJ(1));
	t121 = cos(qJ(2));
	t127 = t119 * t121;
	t124 = qJ(3) + pkin(10);
	t114 = qJ(5) + t124;
	t111 = sin(t114);
	t122 = cos(qJ(1));
	t126 = t122 * t111;
	t112 = cos(t114);
	t125 = t122 * t112;
	t113 = cos(pkin(10)) * pkin(4) + pkin(3);
	t117 = sin(qJ(3));
	t120 = cos(qJ(3));
	t123 = t117 * t113 + t120 * t128 + pkin(7);
	t118 = sin(qJ(2));
	t115 = pkin(9) + qJ(4) + pkin(8);
	t110 = t119 * t111 + t121 * t125;
	t109 = -t119 * t112 + t121 * t126;
	t108 = t112 * t127 - t126;
	t107 = t111 * t127 + t125;
	t106 = (t120 * t113 - t117 * t128 + pkin(2)) * t121 + t115 * t118 + pkin(1);
	t1 = [t110, t122 * t118, t109, t110 * pkin(5) + t109 * qJ(6) + t106 * t122 + t123 * t119 + 0; t108, t119 * t118, t107, t108 * pkin(5) + t107 * qJ(6) + t106 * t119 - t123 * t122 + 0; t118 * t112, -t121, t118 * t111, -t121 * t115 + pkin(6) + 0 + (pkin(4) * cos(t124) + t120 * pkin(3) + t112 * pkin(5) + t111 * qJ(6) + pkin(2)) * t118; 0, 0, 0, 1;];
	Tc_mdh = t1;
end