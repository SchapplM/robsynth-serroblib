% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR14 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for the body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRRR14_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR14_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:42:10
	% EndTime: 2024-09-27 18:42:10
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:42:10
	% EndTime: 2024-09-27 18:42:10
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t66 = cos(qJ(1));
	t65 = sin(qJ(1));
	t1 = [t66, -t65, 0, 0; t65, t66, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:42:10
	% EndTime: 2024-09-27 18:42:10
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t69 = qJ(1) + qJ(2);
	t68 = cos(t69);
	t67 = sin(t69);
	t1 = [t68, -t67, 0, cos(qJ(1)) * pkin(1) + 0; t67, t68, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, pkin(7) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:42:10
	% EndTime: 2024-09-27 18:42:10
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (28->15), mult. (25->21), div. (0->0), fcn. (38->8), ass. (0->12)
	t72 = qJ(1) + qJ(2);
	t70 = sin(t72);
	t73 = sin(pkin(5));
	t80 = t70 * t73;
	t71 = cos(t72);
	t79 = t71 * t73;
	t74 = cos(pkin(5));
	t75 = sin(qJ(3));
	t78 = t74 * t75;
	t76 = cos(qJ(3));
	t77 = t74 * t76;
	t1 = [-t70 * t78 + t71 * t76, -t70 * t77 - t71 * t75, t80, t71 * pkin(2) + pkin(8) * t80 + cos(qJ(1)) * pkin(1) + 0; t70 * t76 + t71 * t78, -t70 * t75 + t71 * t77, -t79, t70 * pkin(2) - pkin(8) * t79 + sin(qJ(1)) * pkin(1) + 0; t73 * t75, t73 * t76, t74, t74 * pkin(8) + pkin(6) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:42:10
	% EndTime: 2024-09-27 18:42:10
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (70->26), mult. (39->28), div. (0->0), fcn. (46->14), ass. (0->21)
	t101 = pkin(3) * sin(qJ(3));
	t95 = qJ(3) + qJ(4);
	t100 = pkin(8) + pkin(9);
	t98 = cos(pkin(5));
	t97 = sin(pkin(5));
	t96 = qJ(1) + qJ(2);
	t94 = cos(t96);
	t93 = cos(t95);
	t92 = sin(t96);
	t91 = sin(t95);
	t90 = pkin(5) - t95;
	t89 = pkin(5) + t95;
	t88 = cos(qJ(3)) * pkin(3) + pkin(2);
	t87 = cos(t89);
	t86 = sin(t90);
	t85 = cos(t90) / 0.2e1;
	t84 = sin(t89) / 0.2e1;
	t83 = -t97 * t100 + t98 * t101;
	t82 = t85 + t87 / 0.2e1;
	t81 = t84 - t86 / 0.2e1;
	t1 = [-t92 * t81 + t94 * t93, -t92 * t82 - t94 * t91, t92 * t97, t94 * t88 - t92 * t83 + cos(qJ(1)) * pkin(1) + 0; t94 * t81 + t92 * t93, t94 * t82 - t92 * t91, -t94 * t97, t92 * t88 + t94 * t83 + sin(qJ(1)) * pkin(1) + 0; t85 - t87 / 0.2e1, t84 + t86 / 0.2e1, t98, t98 * t100 + t97 * t101 + pkin(6) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 18:42:10
	% EndTime: 2024-09-27 18:42:10
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (99->30), mult. (53->32), div. (0->0), fcn. (60->17), ass. (0->23)
	t126 = qJ(3) + qJ(4);
	t117 = qJ(5) + t126;
	t124 = cos(qJ(3));
	t125 = pkin(4) * sin(qJ(4)) * t124 + (cos(qJ(4)) * pkin(4) + pkin(3)) * sin(qJ(3));
	t121 = cos(pkin(5));
	t120 = sin(pkin(5));
	t119 = qJ(1) + qJ(2);
	t118 = pkin(8) + pkin(9) + pkin(10);
	t116 = cos(t119);
	t115 = sin(t119);
	t113 = cos(t117);
	t112 = sin(t117);
	t111 = pkin(5) - t117;
	t110 = pkin(5) + t117;
	t109 = cos(t110);
	t108 = sin(t111);
	t107 = cos(t111) / 0.2e1;
	t106 = sin(t110) / 0.2e1;
	t105 = pkin(2) + pkin(4) * cos(t126) + t124 * pkin(3);
	t104 = t107 + t109 / 0.2e1;
	t103 = t106 - t108 / 0.2e1;
	t102 = -t120 * t118 + t125 * t121;
	t1 = [-t115 * t103 + t116 * t113, -t115 * t104 - t116 * t112, t115 * t120, t116 * t105 - t115 * t102 + cos(qJ(1)) * pkin(1) + 0; t116 * t103 + t115 * t113, t116 * t104 - t115 * t112, -t116 * t120, t115 * t105 + t116 * t102 + sin(qJ(1)) * pkin(1) + 0; t107 - t109 / 0.2e1, t106 + t108 / 0.2e1, t121, t121 * t118 + t125 * t120 + pkin(6) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end