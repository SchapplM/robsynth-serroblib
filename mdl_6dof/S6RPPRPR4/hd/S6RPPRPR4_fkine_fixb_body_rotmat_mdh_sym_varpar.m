% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRPR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:25
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:25:42
	% EndTime: 2020-11-04 21:25:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:25:42
	% EndTime: 2020-11-04 21:25:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t1 = [t56, -t55, 0, 0; t55, t56, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:25:42
	% EndTime: 2020-11-04 21:25:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t58 = cos(qJ(1));
	t57 = sin(qJ(1));
	t1 = [t58, 0, t57, t58 * pkin(1) + t57 * qJ(2) + 0; t57, 0, -t58, t57 * pkin(1) - t58 * qJ(2) + 0; 0, 1, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:25:42
	% EndTime: 2020-11-04 21:25:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t65 = pkin(1) + pkin(2);
	t64 = cos(qJ(1));
	t63 = sin(qJ(1));
	t62 = cos(pkin(9));
	t61 = sin(pkin(9));
	t60 = -t64 * t61 + t63 * t62;
	t59 = -t63 * t61 - t64 * t62;
	t1 = [-t59, t60, 0, t63 * qJ(2) + t65 * t64 + 0; t60, t59, 0, -t64 * qJ(2) + t65 * t63 + 0; 0, 0, -1, -qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:25:42
	% EndTime: 2020-11-04 21:25:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (26->17), mult. (28->16), div. (0->0), fcn. (42->6), ass. (0->11)
	t76 = cos(qJ(1));
	t75 = cos(qJ(4));
	t74 = sin(qJ(1));
	t73 = sin(qJ(4));
	t72 = cos(pkin(9));
	t71 = sin(pkin(9));
	t69 = -t71 * pkin(3) + t72 * pkin(7) - qJ(2);
	t68 = t72 * pkin(3) + t71 * pkin(7) + pkin(1) + pkin(2);
	t67 = t74 * t71 + t76 * t72;
	t66 = t76 * t71 - t74 * t72;
	t1 = [t67 * t75, -t67 * t73, t66, t68 * t76 - t69 * t74 + 0; -t66 * t75, t66 * t73, t67, t68 * t74 + t69 * t76 + 0; -t73, -t75, 0, -qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:25:42
	% EndTime: 2020-11-04 21:25:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (41->21), mult. (37->18), div. (0->0), fcn. (51->8), ass. (0->14)
	t84 = sin(pkin(9));
	t85 = cos(pkin(9));
	t86 = qJ(5) + pkin(7);
	t91 = pkin(4) * cos(qJ(4)) + pkin(3);
	t92 = t91 * t84 - t86 * t85 + qJ(2);
	t89 = cos(qJ(1));
	t87 = sin(qJ(1));
	t83 = qJ(4) + pkin(10);
	t82 = cos(t83);
	t81 = sin(t83);
	t79 = t87 * t84 + t89 * t85;
	t78 = t89 * t84 - t87 * t85;
	t77 = t86 * t84 + t91 * t85 + pkin(1) + pkin(2);
	t1 = [t79 * t82, -t79 * t81, t78, t77 * t89 + t92 * t87 + 0; -t78 * t82, t78 * t81, t79, t77 * t87 - t92 * t89 + 0; -t81, -t82, 0, -sin(qJ(4)) * pkin(4) - qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:25:42
	% EndTime: 2020-11-04 21:25:42
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (65->28), mult. (75->32), div. (0->0), fcn. (102->10), ass. (0->21)
	t100 = sin(pkin(9));
	t101 = cos(pkin(9));
	t102 = qJ(5) + pkin(7);
	t109 = pkin(4) * cos(qJ(4)) + pkin(3);
	t115 = t109 * t100 - t102 * t101 + qJ(2);
	t114 = sin(qJ(1));
	t103 = sin(qJ(6));
	t106 = cos(qJ(1));
	t94 = t106 * t100 - t114 * t101;
	t113 = t94 * t103;
	t104 = cos(qJ(6));
	t112 = t94 * t104;
	t95 = t114 * t100 + t106 * t101;
	t111 = t95 * t103;
	t110 = t95 * t104;
	t99 = qJ(4) + pkin(10);
	t97 = sin(t99);
	t98 = cos(t99);
	t107 = t98 * pkin(5) + t97 * pkin(8);
	t93 = t102 * t100 + t109 * t101 + pkin(1) + pkin(2);
	t1 = [t98 * t110 + t113, -t98 * t111 + t112, t95 * t97, t93 * t106 + t107 * t95 + t115 * t114 + 0; -t98 * t112 + t111, t98 * t113 + t110, -t94 * t97, -t115 * t106 - t107 * t94 + t93 * t114 + 0; -t97 * t104, t97 * t103, t98, -t97 * pkin(5) + t98 * pkin(8) - sin(qJ(4)) * pkin(4) - qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end