% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRRP3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:54
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:28
	% EndTime: 2020-11-04 19:54:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:28
	% EndTime: 2020-11-04 19:54:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t66 = cos(pkin(7));
	t65 = sin(pkin(7));
	t1 = [t66, -t65, 0, 0; t65, t66, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:28
	% EndTime: 2020-11-04 19:54:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t70 = cos(pkin(7));
	t69 = cos(pkin(8));
	t68 = sin(pkin(7));
	t67 = sin(pkin(8));
	t1 = [t70 * t69, -t70 * t67, t68, t70 * pkin(1) + t68 * qJ(2) + 0; t68 * t69, -t68 * t67, -t70, t68 * pkin(1) - t70 * qJ(2) + 0; t67, t69, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:28
	% EndTime: 2020-11-04 19:54:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->12)
	t73 = sin(pkin(7));
	t76 = sin(qJ(3));
	t81 = t73 * t76;
	t77 = cos(qJ(3));
	t80 = t73 * t77;
	t75 = cos(pkin(7));
	t79 = t75 * t76;
	t78 = t75 * t77;
	t74 = cos(pkin(8));
	t72 = sin(pkin(8));
	t71 = t74 * pkin(2) + t72 * pkin(5) + pkin(1);
	t1 = [t74 * t78 + t81, -t74 * t79 + t80, t75 * t72, t73 * qJ(2) + t71 * t75 + 0; t74 * t80 - t79, -t74 * t81 - t78, t73 * t72, -t75 * qJ(2) + t71 * t73 + 0; t72 * t77, -t72 * t76, -t74, t72 * pkin(2) - t74 * pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:28
	% EndTime: 2020-11-04 19:54:28
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (33->29), mult. (67->48), div. (0->0), fcn. (88->8), ass. (0->21)
	t85 = sin(pkin(8));
	t89 = sin(qJ(4));
	t101 = t85 * t89;
	t91 = cos(qJ(4));
	t100 = t85 * t91;
	t92 = cos(qJ(3));
	t99 = t85 * t92;
	t86 = sin(pkin(7));
	t87 = cos(pkin(8));
	t98 = t86 * t87;
	t90 = sin(qJ(3));
	t97 = t86 * t90;
	t96 = t86 * t92;
	t88 = cos(pkin(7));
	t95 = t87 * t88;
	t94 = t88 * t90;
	t93 = t88 * t92;
	t84 = t87 * pkin(2) + t85 * pkin(5) + pkin(1);
	t83 = t87 * t93 + t97;
	t82 = t87 * t96 - t94;
	t1 = [t88 * t101 + t83 * t91, t88 * t100 - t83 * t89, t87 * t94 - t96, (pkin(3) * t95 - t86 * pkin(6)) * t92 + (t86 * pkin(3) + pkin(6) * t95) * t90 + t84 * t88 + t86 * qJ(2) + 0; t86 * t101 + t82 * t91, t86 * t100 - t82 * t89, t87 * t97 + t93, (pkin(3) * t98 + t88 * pkin(6)) * t92 + (-t88 * pkin(3) + pkin(6) * t98) * t90 + t84 * t86 - t88 * qJ(2) + 0; -t87 * t89 + t91 * t99, -t87 * t91 - t89 * t99, t85 * t90, -t87 * pkin(5) + qJ(1) + 0 + (pkin(3) * t92 + pkin(6) * t90 + pkin(2)) * t85; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:28
	% EndTime: 2020-11-04 19:54:28
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (49->27), mult. (105->39), div. (0->0), fcn. (126->8), ass. (0->22)
	t110 = sin(qJ(3));
	t112 = cos(qJ(3));
	t109 = sin(qJ(4));
	t111 = cos(qJ(4));
	t115 = t111 * pkin(4) + qJ(5) * t109 + pkin(3);
	t126 = pkin(6) * t110 + t115 * t112 + pkin(2);
	t125 = t109 * pkin(4) - qJ(5) * t111 + pkin(5);
	t124 = -pkin(6) * t112 + t115 * t110 + qJ(2);
	t106 = sin(pkin(7));
	t123 = t106 * t110;
	t122 = t106 * t112;
	t108 = cos(pkin(7));
	t121 = t108 * t110;
	t120 = t108 * t112;
	t105 = sin(pkin(8));
	t119 = t109 * t105;
	t118 = t111 * t105;
	t107 = cos(pkin(8));
	t113 = t125 * t105 + t126 * t107 + pkin(1);
	t103 = t107 * t120 + t123;
	t102 = t107 * t122 - t121;
	t1 = [t103 * t111 + t108 * t119, t107 * t121 - t122, t103 * t109 - t108 * t118, t124 * t106 + t113 * t108 + 0; t102 * t111 + t106 * t119, t107 * t123 + t120, t102 * t109 - t106 * t118, t113 * t106 - t124 * t108 + 0; -t107 * t109 + t112 * t118, t105 * t110, t107 * t111 + t112 * t119, t126 * t105 - t125 * t107 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end