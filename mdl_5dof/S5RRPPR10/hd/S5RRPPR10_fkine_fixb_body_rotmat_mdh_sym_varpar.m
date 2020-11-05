% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPR10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:31
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPPR10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:57
	% EndTime: 2020-11-04 20:31:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:57
	% EndTime: 2020-11-04 20:31:57
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t49 = cos(qJ(1));
	t48 = sin(qJ(1));
	t1 = [t49, -t48, 0, 0; t48, t49, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:57
	% EndTime: 2020-11-04 20:31:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t53 = cos(qJ(1));
	t52 = cos(qJ(2));
	t51 = sin(qJ(1));
	t50 = sin(qJ(2));
	t1 = [t53 * t52, -t53 * t50, t51, t53 * pkin(1) + t51 * pkin(6) + 0; t51 * t52, -t51 * t50, -t53, t51 * pkin(1) - t53 * pkin(6) + 0; t50, t52, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:57
	% EndTime: 2020-11-04 20:31:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t58 = sin(qJ(1));
	t59 = cos(qJ(2));
	t63 = t58 * t59;
	t55 = sin(pkin(8));
	t60 = cos(qJ(1));
	t62 = t60 * t55;
	t56 = cos(pkin(8));
	t61 = t60 * t56;
	t57 = sin(qJ(2));
	t54 = t59 * pkin(2) + t57 * qJ(3) + pkin(1);
	t1 = [t58 * t55 + t59 * t61, t58 * t56 - t59 * t62, t60 * t57, t58 * pkin(6) + t54 * t60 + 0; t56 * t63 - t62, -t55 * t63 - t61, t58 * t57, -t60 * pkin(6) + t54 * t58 + 0; t57 * t56, -t57 * t55, -t59, t57 * pkin(2) - t59 * qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:57
	% EndTime: 2020-11-04 20:31:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (26->18), mult. (36->25), div. (0->0), fcn. (49->6), ass. (0->13)
	t70 = sin(qJ(1));
	t71 = cos(qJ(2));
	t75 = t70 * t71;
	t67 = sin(pkin(8));
	t72 = cos(qJ(1));
	t74 = t72 * t67;
	t68 = cos(pkin(8));
	t73 = t72 * t68;
	t69 = sin(qJ(2));
	t66 = -t67 * pkin(3) + qJ(4) * t68 - pkin(6);
	t65 = t68 * pkin(3) + t67 * qJ(4) + pkin(2);
	t64 = t69 * qJ(3) + t65 * t71 + pkin(1);
	t1 = [t70 * t67 + t71 * t73, t72 * t69, -t70 * t68 + t71 * t74, t64 * t72 - t66 * t70 + 0; t68 * t75 - t74, t70 * t69, t67 * t75 + t73, t64 * t70 + t66 * t72 + 0; t69 * t68, -t71, t69 * t67, -t71 * qJ(3) + t65 * t69 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:31:57
	% EndTime: 2020-11-04 20:31:57
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (45->23), mult. (56->30), div. (0->0), fcn. (79->8), ass. (0->18)
	t85 = sin(qJ(1));
	t87 = cos(qJ(2));
	t92 = t85 * t87;
	t88 = cos(qJ(1));
	t91 = t87 * t88;
	t80 = sin(pkin(8));
	t81 = cos(pkin(8));
	t89 = pkin(3) + pkin(4);
	t90 = qJ(4) * t81 - t89 * t80 - pkin(6);
	t86 = cos(qJ(5));
	t84 = sin(qJ(2));
	t83 = sin(qJ(5));
	t82 = qJ(3) - pkin(7);
	t79 = t80 * qJ(4) + t89 * t81 + pkin(2);
	t78 = t83 * t80 + t86 * t81;
	t77 = t86 * t80 - t83 * t81;
	t76 = t79 * t87 + t82 * t84 + pkin(1);
	t1 = [t85 * t77 + t78 * t91, t77 * t91 - t85 * t78, -t88 * t84, t76 * t88 - t90 * t85 + 0; -t77 * t88 + t78 * t92, t77 * t92 + t78 * t88, -t85 * t84, t76 * t85 + t90 * t88 + 0; t84 * t78, t84 * t77, t87, t79 * t84 - t82 * t87 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end