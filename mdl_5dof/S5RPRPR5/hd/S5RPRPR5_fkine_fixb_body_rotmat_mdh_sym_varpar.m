% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:19
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:19:21
	% EndTime: 2020-11-04 20:19:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:19:21
	% EndTime: 2020-11-04 20:19:21
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t51 = cos(qJ(1));
	t50 = sin(qJ(1));
	t1 = [0, 0, 1, pkin(5) + 0; t50, t51, 0, 0; -t51, t50, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:19:21
	% EndTime: 2020-11-04 20:19:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->9), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t55 = cos(qJ(1));
	t54 = sin(qJ(1));
	t53 = cos(pkin(8));
	t52 = sin(pkin(8));
	t1 = [t52, t53, 0, pkin(5) + 0; t54 * t53, -t54 * t52, -t55, t54 * pkin(1) - t55 * qJ(2) + 0; -t55 * t53, t55 * t52, -t54, -t55 * pkin(1) - t54 * qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:19:21
	% EndTime: 2020-11-04 20:19:21
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (18->16), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->12)
	t59 = sin(qJ(3));
	t60 = sin(qJ(1));
	t66 = t60 * t59;
	t61 = cos(qJ(3));
	t65 = t60 * t61;
	t62 = cos(qJ(1));
	t64 = t62 * t59;
	t63 = t62 * t61;
	t58 = cos(pkin(8));
	t57 = sin(pkin(8));
	t56 = pkin(2) * t58 + t57 * pkin(6) + pkin(1);
	t1 = [t57 * t61, -t57 * t59, -t58, t57 * pkin(2) - t58 * pkin(6) + pkin(5) + 0; t58 * t65 - t64, -t58 * t66 - t63, t60 * t57, -t62 * qJ(2) + t56 * t60 + 0; -t58 * t63 - t66, t58 * t64 - t65, -t62 * t57, -t60 * qJ(2) - t56 * t62 + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:19:21
	% EndTime: 2020-11-04 20:19:21
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (36->20), mult. (34->22), div. (0->0), fcn. (47->8), ass. (0->16)
	t71 = sin(pkin(8));
	t72 = cos(pkin(8));
	t73 = qJ(4) + pkin(6);
	t78 = pkin(3) * cos(qJ(3)) + pkin(2);
	t83 = t73 * t71 + t78 * t72 + pkin(1);
	t70 = qJ(3) + pkin(9);
	t68 = sin(t70);
	t74 = sin(qJ(1));
	t82 = t74 * t68;
	t69 = cos(t70);
	t81 = t74 * t69;
	t76 = cos(qJ(1));
	t80 = t76 * t68;
	t79 = t76 * t69;
	t67 = -sin(qJ(3)) * pkin(3) - qJ(2);
	t1 = [t71 * t69, -t71 * t68, -t72, t78 * t71 - t72 * t73 + pkin(5) + 0; t72 * t81 - t80, -t72 * t82 - t79, t74 * t71, t67 * t76 + t83 * t74 + 0; -t72 * t79 - t82, t72 * t80 - t81, -t76 * t71, t67 * t74 - t83 * t76 + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:19:21
	% EndTime: 2020-11-04 20:19:21
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (59->23), mult. (42->24), div. (0->0), fcn. (55->10), ass. (0->17)
	t90 = qJ(3) + pkin(9);
	t88 = qJ(5) + t90;
	t86 = sin(t88);
	t93 = sin(qJ(1));
	t100 = t93 * t86;
	t87 = cos(t88);
	t99 = t93 * t87;
	t94 = cos(qJ(1));
	t98 = t94 * t86;
	t97 = t94 * t87;
	t96 = -qJ(2) - pkin(4) * sin(t90) - sin(qJ(3)) * pkin(3);
	t84 = pkin(4) * cos(t90) + pkin(3) * cos(qJ(3)) + pkin(2);
	t89 = -pkin(7) - qJ(4) - pkin(6);
	t91 = sin(pkin(8));
	t92 = cos(pkin(8));
	t95 = t84 * t92 - t89 * t91 + pkin(1);
	t1 = [t91 * t87, -t91 * t86, -t92, t91 * t84 + t92 * t89 + pkin(5) + 0; t92 * t99 - t98, -t92 * t100 - t97, t93 * t91, t95 * t93 + t96 * t94 + 0; -t92 * t97 - t100, t92 * t98 - t99, -t94 * t91, t96 * t93 - t95 * t94 + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end