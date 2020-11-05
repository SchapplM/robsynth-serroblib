% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPPRR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:52
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:52:15
	% EndTime: 2020-11-04 19:52:15
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:52:15
	% EndTime: 2020-11-04 19:52:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t64 = cos(pkin(7));
	t63 = sin(pkin(7));
	t1 = [t64, -t63, 0, 0; t63, t64, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:52:15
	% EndTime: 2020-11-04 19:52:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t68 = cos(pkin(7));
	t67 = cos(pkin(8));
	t66 = sin(pkin(7));
	t65 = sin(pkin(8));
	t1 = [t68 * t67, -t68 * t65, t66, t68 * pkin(1) + t66 * qJ(2) + 0; t66 * t67, -t66 * t65, -t68, t66 * pkin(1) - t68 * qJ(2) + 0; t65, t67, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:52:15
	% EndTime: 2020-11-04 19:52:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t72 = sin(pkin(7));
	t74 = cos(pkin(8));
	t78 = t72 * t74;
	t70 = sin(pkin(9));
	t75 = cos(pkin(7));
	t77 = t75 * t70;
	t73 = cos(pkin(9));
	t76 = t75 * t73;
	t71 = sin(pkin(8));
	t69 = t74 * pkin(2) + t71 * qJ(3) + pkin(1);
	t1 = [t72 * t70 + t74 * t76, t72 * t73 - t74 * t77, t75 * t71, t72 * qJ(2) + t69 * t75 + 0; t73 * t78 - t77, -t70 * t78 - t76, t72 * t71, -t75 * qJ(2) + t69 * t72 + 0; t71 * t73, -t71 * t70, -t74, t71 * pkin(2) - t74 * qJ(3) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:52:15
	% EndTime: 2020-11-04 19:52:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (33->23), mult. (57->36), div. (0->0), fcn. (78->8), ass. (0->19)
	t85 = sin(pkin(8));
	t90 = sin(qJ(4));
	t96 = t85 * t90;
	t91 = cos(qJ(4));
	t95 = t85 * t91;
	t86 = sin(pkin(7));
	t88 = cos(pkin(8));
	t94 = t86 * t88;
	t84 = sin(pkin(9));
	t89 = cos(pkin(7));
	t93 = t89 * t84;
	t87 = cos(pkin(9));
	t92 = t89 * t87;
	t83 = t87 * pkin(3) + t84 * pkin(5) + pkin(2);
	t82 = -t84 * pkin(3) + t87 * pkin(5) - qJ(2);
	t81 = t86 * t84 + t88 * t92;
	t80 = t87 * t94 - t93;
	t79 = t85 * qJ(3) + t83 * t88 + pkin(1);
	t1 = [t81 * t91 + t89 * t96, -t81 * t90 + t89 * t95, -t86 * t87 + t88 * t93, t79 * t89 - t82 * t86 + 0; t80 * t91 + t86 * t96, -t80 * t90 + t86 * t95, t84 * t94 + t92, t79 * t86 + t82 * t89 + 0; t87 * t95 - t88 * t90, -t87 * t96 - t88 * t91, t85 * t84, -t88 * qJ(3) + t83 * t85 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:52:15
	% EndTime: 2020-11-04 19:52:16
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (62->43), mult. (127->71), div. (0->0), fcn. (161->10), ass. (0->31)
	t108 = sin(pkin(9));
	t109 = sin(pkin(8));
	t127 = t108 * t109;
	t111 = cos(pkin(9));
	t126 = t109 * t111;
	t115 = sin(qJ(4));
	t125 = t109 * t115;
	t117 = cos(qJ(4));
	t124 = t109 * t117;
	t110 = sin(pkin(7));
	t123 = t110 * t108;
	t112 = cos(pkin(8));
	t122 = t111 * t112;
	t113 = cos(pkin(7));
	t121 = t113 * t108;
	t120 = t113 * t111;
	t99 = t110 * t122 - t121;
	t119 = t110 * t125 + t99 * t117;
	t101 = t112 * t120 + t123;
	t118 = t101 * t117 + t113 * t125;
	t116 = cos(qJ(5));
	t114 = sin(qJ(5));
	t106 = t111 * pkin(3) + t108 * pkin(5) + pkin(2);
	t105 = -t108 * pkin(3) + t111 * pkin(5) - qJ(2);
	t104 = pkin(4) * t122 - t109 * pkin(6);
	t103 = t109 * pkin(4) + pkin(6) * t122;
	t102 = t111 * t124 - t112 * t115;
	t100 = t110 * t111 - t112 * t121;
	t98 = t112 * t123 + t120;
	t97 = t109 * qJ(3) + t106 * t112 + pkin(1);
	t1 = [-t114 * t100 + t118 * t116, -t116 * t100 - t118 * t114, t101 * t115 - t113 * t124, (pkin(4) * t123 + t104 * t113) * t117 + (pkin(6) * t123 + t103 * t113) * t115 + t97 * t113 - t105 * t110 + 0; t98 * t114 + t119 * t116, -t119 * t114 + t98 * t116, -t110 * t124 + t99 * t115, (-pkin(4) * t121 + t104 * t110) * t117 + (-pkin(6) * t121 + t103 * t110) * t115 + t97 * t110 + t105 * t113 + 0; t102 * t116 + t114 * t127, -t102 * t114 + t116 * t127, t111 * t125 + t112 * t117, (pkin(4) * t126 + t112 * pkin(6)) * t117 + (-t112 * pkin(4) + pkin(6) * t126) * t115 + t106 * t109 - t112 * qJ(3) + 0 + qJ(1); 0, 0, 0, 1;];
	Tc_mdh = t1;
end