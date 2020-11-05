% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPP1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:29
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPPP1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:29:13
	% EndTime: 2020-11-04 20:29:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:29:13
	% EndTime: 2020-11-04 20:29:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t1 = [t56, -t55, 0, 0; t55, t56, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:29:13
	% EndTime: 2020-11-04 20:29:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t60 = cos(qJ(1));
	t59 = cos(qJ(2));
	t58 = sin(qJ(1));
	t57 = sin(qJ(2));
	t1 = [t60 * t59, -t60 * t57, t58, pkin(1) * t60 + pkin(7) * t58 + 0; t58 * t59, -t58 * t57, -t60, pkin(1) * t58 - pkin(7) * t60 + 0; t57, t59, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:29:13
	% EndTime: 2020-11-04 20:29:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (26->21), mult. (56->34), div. (0->0), fcn. (77->8), ass. (0->21)
	t64 = sin(pkin(5));
	t68 = sin(qJ(1));
	t80 = t64 * t68;
	t70 = cos(qJ(1));
	t79 = t64 * t70;
	t66 = cos(pkin(5));
	t69 = cos(qJ(2));
	t78 = t66 * t69;
	t77 = t68 * t66;
	t76 = t68 * t69;
	t75 = t69 * t70;
	t74 = t70 * t66;
	t73 = qJ(3) * t64;
	t67 = sin(qJ(2));
	t72 = -t67 * t74 + t80;
	t71 = -t67 * t77 - t79;
	t65 = cos(pkin(8));
	t63 = sin(pkin(8));
	t62 = t66 * qJ(3) + pkin(7);
	t61 = pkin(2) * t69 + t67 * t73 + pkin(1);
	t1 = [t72 * t63 + t65 * t75, -t63 * t75 + t72 * t65, t67 * t79 + t77, t61 * t70 + t62 * t68 + 0; t71 * t63 + t65 * t76, -t63 * t76 + t71 * t65, t67 * t80 - t74, t61 * t68 - t62 * t70 + 0; t63 * t78 + t67 * t65, -t67 * t63 + t65 * t78, -t69 * t64, t67 * pkin(2) - t69 * t73 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:29:13
	% EndTime: 2020-11-04 20:29:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (42->26), mult. (79->40), div. (0->0), fcn. (100->8), ass. (0->23)
	t86 = sin(pkin(5));
	t90 = sin(qJ(1));
	t102 = t86 * t90;
	t92 = cos(qJ(1));
	t101 = t86 * t92;
	t88 = cos(pkin(5));
	t91 = cos(qJ(2));
	t100 = t88 * t91;
	t99 = t90 * t88;
	t98 = t90 * t91;
	t97 = t91 * t92;
	t96 = t92 * t88;
	t85 = sin(pkin(8));
	t87 = cos(pkin(8));
	t84 = -t85 * pkin(3) + t87 * qJ(4);
	t95 = t86 * qJ(3) + t84 * t88;
	t89 = sin(qJ(2));
	t94 = t89 * t96 - t102;
	t93 = t89 * t99 + t101;
	t83 = t87 * pkin(3) + t85 * qJ(4) + pkin(2);
	t82 = -t88 * qJ(3) + t84 * t86 - pkin(7);
	t81 = t83 * t91 + t95 * t89 + pkin(1);
	t1 = [t89 * t101 + t99, t94 * t85 - t87 * t97, t85 * t97 + t94 * t87, t81 * t92 - t82 * t90 + 0; t89 * t102 - t96, t93 * t85 - t87 * t98, t85 * t98 + t93 * t87, t81 * t90 + t82 * t92 + 0; -t91 * t86, -t85 * t100 - t89 * t87, -t87 * t100 + t89 * t85, t83 * t89 - t95 * t91 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:29:13
	% EndTime: 2020-11-04 20:29:14
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (55->28), mult. (77->40), div. (0->0), fcn. (98->8), ass. (0->25)
	t107 = sin(pkin(8));
	t109 = cos(pkin(8));
	t111 = qJ(5) + pkin(3);
	t106 = t109 * qJ(4) - t107 * t111;
	t108 = sin(pkin(5));
	t110 = cos(pkin(5));
	t112 = -qJ(3) - pkin(4);
	t129 = t106 * t108 + t112 * t110 - pkin(7);
	t128 = t106 * t110 - t108 * t112;
	t114 = sin(qJ(1));
	t126 = t108 * t114;
	t116 = cos(qJ(1));
	t125 = t108 * t116;
	t115 = cos(qJ(2));
	t124 = t110 * t115;
	t123 = t114 * t110;
	t122 = t114 * t115;
	t121 = t115 * t116;
	t120 = t116 * t110;
	t113 = sin(qJ(2));
	t118 = t113 * t120 - t126;
	t117 = t113 * t123 + t125;
	t104 = t107 * qJ(4) + t111 * t109 + pkin(2);
	t103 = t104 * t115 + t128 * t113 + pkin(1);
	t1 = [t113 * t125 + t123, t107 * t121 + t118 * t109, -t118 * t107 + t109 * t121, t103 * t116 - t129 * t114 + 0; t113 * t126 - t120, t107 * t122 + t117 * t109, -t117 * t107 + t109 * t122, t103 * t114 + t129 * t116 + 0; -t115 * t108, t113 * t107 - t109 * t124, t107 * t124 + t113 * t109, t104 * t113 - t128 * t115 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end