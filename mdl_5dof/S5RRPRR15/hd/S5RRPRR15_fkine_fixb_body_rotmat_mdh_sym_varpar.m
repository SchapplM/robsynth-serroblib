% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR15 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:39
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPRR15_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR15_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:39:05
	% EndTime: 2020-11-04 20:39:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:39:05
	% EndTime: 2020-11-04 20:39:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t40 = cos(qJ(1));
	t39 = sin(qJ(1));
	t1 = [t40, -t39, 0, 0; t39, t40, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:39:05
	% EndTime: 2020-11-04 20:39:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t44 = cos(qJ(1));
	t43 = cos(qJ(2));
	t42 = sin(qJ(1));
	t41 = sin(qJ(2));
	t1 = [t44 * t43, -t44 * t41, t42, t44 * pkin(1) + t42 * pkin(6) + 0; t42 * t43, -t42 * t41, -t44, t42 * pkin(1) - t44 * pkin(6) + 0; t41, t43, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:39:05
	% EndTime: 2020-11-04 20:39:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (16->14), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t49 = cos(qJ(1));
	t48 = cos(qJ(2));
	t47 = sin(qJ(1));
	t46 = sin(qJ(2));
	t45 = t48 * pkin(2) + t46 * qJ(3) + pkin(1);
	t1 = [t47, -t49 * t48, t49 * t46, t47 * pkin(6) + t45 * t49 + 0; -t49, -t47 * t48, t47 * t46, -t49 * pkin(6) + t45 * t47 + 0; 0, -t46, -t48, t46 * pkin(2) - t48 * qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:39:05
	% EndTime: 2020-11-04 20:39:05
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->14)
	t51 = sin(qJ(4));
	t53 = sin(qJ(1));
	t62 = t53 * t51;
	t54 = cos(qJ(4));
	t61 = t53 * t54;
	t56 = cos(qJ(1));
	t60 = t56 * t51;
	t59 = t56 * t54;
	t58 = pkin(2) + pkin(7);
	t57 = pkin(3) + pkin(6);
	t55 = cos(qJ(2));
	t52 = sin(qJ(2));
	t50 = t52 * qJ(3) + t58 * t55 + pkin(1);
	t1 = [t52 * t60 + t61, t52 * t59 - t62, t56 * t55, t50 * t56 + t57 * t53 + 0; t52 * t62 - t59, t52 * t61 + t60, t53 * t55, t50 * t53 - t57 * t56 + 0; -t55 * t51, -t55 * t54, t52, -t55 * qJ(3) + t58 * t52 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:39:05
	% EndTime: 2020-11-04 20:39:05
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (40->20), mult. (31->22), div. (0->0), fcn. (44->8), ass. (0->16)
	t69 = qJ(4) + qJ(5);
	t66 = sin(t69);
	t71 = sin(qJ(1));
	t77 = t71 * t66;
	t67 = cos(t69);
	t76 = t71 * t67;
	t73 = cos(qJ(1));
	t75 = t73 * t66;
	t74 = t73 * t67;
	t72 = cos(qJ(2));
	t70 = sin(qJ(2));
	t68 = pkin(2) + pkin(7) + pkin(8);
	t65 = sin(qJ(4)) * pkin(4) + qJ(3);
	t64 = cos(qJ(4)) * pkin(4) + pkin(3) + pkin(6);
	t63 = t65 * t70 + t68 * t72 + pkin(1);
	t1 = [t70 * t75 + t76, t70 * t74 - t77, t73 * t72, t63 * t73 + t64 * t71 + 0; t70 * t77 - t74, t70 * t76 + t75, t71 * t72, t63 * t71 - t64 * t73 + 0; -t72 * t66, -t72 * t67, t70, -t65 * t72 + t68 * t70 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end