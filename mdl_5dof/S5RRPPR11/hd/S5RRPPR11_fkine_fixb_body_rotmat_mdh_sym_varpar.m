% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPR11 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:32
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPPR11_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR11_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:32:14
	% EndTime: 2020-11-04 20:32:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:32:14
	% EndTime: 2020-11-04 20:32:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t38 = cos(qJ(1));
	t37 = sin(qJ(1));
	t1 = [t38, -t37, 0, 0; t37, t38, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:32:14
	% EndTime: 2020-11-04 20:32:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t42 = cos(qJ(1));
	t41 = cos(qJ(2));
	t40 = sin(qJ(1));
	t39 = sin(qJ(2));
	t1 = [t42 * t41, -t42 * t39, t40, t42 * pkin(1) + t40 * pkin(6) + 0; t40 * t41, -t40 * t39, -t42, t40 * pkin(1) - t42 * pkin(6) + 0; t39, t41, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:32:14
	% EndTime: 2020-11-04 20:32:14
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (16->14), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t47 = cos(qJ(1));
	t46 = cos(qJ(2));
	t45 = sin(qJ(1));
	t44 = sin(qJ(2));
	t43 = t46 * pkin(2) + t44 * qJ(3) + pkin(1);
	t1 = [t45, -t47 * t46, t47 * t44, t45 * pkin(6) + t43 * t47 + 0; -t47, -t45 * t46, t45 * t44, -t47 * pkin(6) + t43 * t45 + 0; 0, -t44, -t46, t44 * pkin(2) - t46 * qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:32:14
	% EndTime: 2020-11-04 20:32:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->14)
	t49 = sin(pkin(8));
	t53 = sin(qJ(1));
	t60 = t53 * t49;
	t50 = cos(pkin(8));
	t59 = t53 * t50;
	t55 = cos(qJ(1));
	t58 = t55 * t49;
	t57 = t55 * t50;
	t56 = pkin(3) + pkin(6);
	t54 = cos(qJ(2));
	t52 = sin(qJ(2));
	t51 = pkin(2) + qJ(4);
	t48 = t52 * qJ(3) + t51 * t54 + pkin(1);
	t1 = [t52 * t58 + t59, t52 * t57 - t60, t55 * t54, t48 * t55 + t56 * t53 + 0; t52 * t60 - t57, t52 * t59 + t58, t53 * t54, t48 * t53 - t56 * t55 + 0; -t54 * t49, -t54 * t50, t52, -t54 * qJ(3) + t51 * t52 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:32:14
	% EndTime: 2020-11-04 20:32:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (40->20), mult. (31->22), div. (0->0), fcn. (44->8), ass. (0->16)
	t67 = pkin(8) + qJ(5);
	t64 = sin(t67);
	t69 = sin(qJ(1));
	t75 = t69 * t64;
	t65 = cos(t67);
	t74 = t69 * t65;
	t71 = cos(qJ(1));
	t73 = t71 * t64;
	t72 = t71 * t65;
	t70 = cos(qJ(2));
	t68 = sin(qJ(2));
	t66 = pkin(2) + pkin(7) + qJ(4);
	t63 = sin(pkin(8)) * pkin(4) + qJ(3);
	t62 = cos(pkin(8)) * pkin(4) + pkin(3) + pkin(6);
	t61 = t63 * t68 + t66 * t70 + pkin(1);
	t1 = [t68 * t73 + t74, t68 * t72 - t75, t71 * t70, t61 * t71 + t62 * t69 + 0; t68 * t75 - t72, t68 * t74 + t73, t69 * t70, t61 * t69 - t62 * t71 + 0; -t70 * t64, -t70 * t65, t68, -t63 * t70 + t66 * t68 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end