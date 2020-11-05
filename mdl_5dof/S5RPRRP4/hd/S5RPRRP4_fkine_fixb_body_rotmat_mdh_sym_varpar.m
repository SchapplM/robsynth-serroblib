% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:23
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:25
	% EndTime: 2020-11-04 20:23:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:25
	% EndTime: 2020-11-04 20:23:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t52 = cos(qJ(1));
	t51 = sin(qJ(1));
	t1 = [0, 0, 1, pkin(5) + 0; t51, t52, 0, 0; -t52, t51, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:25
	% EndTime: 2020-11-04 20:23:25
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (9->9), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t54 = cos(pkin(8));
	t53 = sin(pkin(8));
	t1 = [t53, t54, 0, pkin(5) + 0; t55 * t54, -t55 * t53, -t56, t55 * pkin(1) - t56 * qJ(2) + 0; -t56 * t54, t56 * t53, -t55, -t56 * pkin(1) - t55 * qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:25
	% EndTime: 2020-11-04 20:23:25
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (18->16), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->12)
	t60 = sin(qJ(3));
	t61 = sin(qJ(1));
	t67 = t61 * t60;
	t62 = cos(qJ(3));
	t66 = t61 * t62;
	t63 = cos(qJ(1));
	t65 = t63 * t60;
	t64 = t63 * t62;
	t59 = cos(pkin(8));
	t58 = sin(pkin(8));
	t57 = pkin(2) * t59 + t58 * pkin(6) + pkin(1);
	t1 = [t58 * t62, -t58 * t60, -t59, t58 * pkin(2) - t59 * pkin(6) + pkin(5) + 0; t59 * t66 - t65, -t59 * t67 - t64, t61 * t58, -t63 * qJ(2) + t57 * t61 + 0; -t59 * t64 - t67, t59 * t65 - t66, -t63 * t58, -t61 * qJ(2) - t57 * t63 + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:25
	% EndTime: 2020-11-04 20:23:25
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (36->20), mult. (34->22), div. (0->0), fcn. (47->8), ass. (0->16)
	t72 = sin(pkin(8));
	t73 = cos(pkin(8));
	t77 = pkin(7) + pkin(6);
	t79 = pkin(3) * cos(qJ(3)) + pkin(2);
	t84 = t77 * t72 + t79 * t73 + pkin(1);
	t71 = qJ(3) + qJ(4);
	t69 = sin(t71);
	t74 = sin(qJ(1));
	t83 = t74 * t69;
	t70 = cos(t71);
	t82 = t74 * t70;
	t76 = cos(qJ(1));
	t81 = t76 * t69;
	t80 = t76 * t70;
	t68 = -sin(qJ(3)) * pkin(3) - qJ(2);
	t1 = [t72 * t70, -t72 * t69, -t73, t79 * t72 - t73 * t77 + pkin(5) + 0; t73 * t82 - t81, -t73 * t83 - t80, t74 * t72, t68 * t76 + t84 * t74 + 0; -t73 * t80 - t83, t73 * t81 - t82, -t76 * t72, t68 * t74 - t84 * t76 + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:25
	% EndTime: 2020-11-04 20:23:25
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (49->22), mult. (42->24), div. (0->0), fcn. (55->8), ass. (0->16)
	t90 = qJ(3) + qJ(4);
	t87 = sin(t90);
	t93 = sin(qJ(1));
	t100 = t93 * t87;
	t88 = cos(t90);
	t99 = t93 * t88;
	t94 = cos(qJ(1));
	t98 = t94 * t87;
	t97 = t94 * t88;
	t96 = -qJ(2) - pkin(4) * t87 - sin(qJ(3)) * pkin(3);
	t85 = pkin(4) * t88 + pkin(3) * cos(qJ(3)) + pkin(2);
	t89 = -qJ(5) - pkin(7) - pkin(6);
	t91 = sin(pkin(8));
	t92 = cos(pkin(8));
	t95 = t85 * t92 - t89 * t91 + pkin(1);
	t1 = [t91 * t88, -t91 * t87, -t92, t91 * t85 + t92 * t89 + pkin(5) + 0; t92 * t99 - t98, -t92 * t100 - t97, t93 * t91, t95 * t93 + t96 * t94 + 0; -t92 * t97 - t100, t92 * t98 - t99, -t94 * t91, t96 * t93 - t95 * t94 + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end