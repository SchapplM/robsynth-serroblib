% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP12 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:25
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRRP12_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP12_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:36
	% EndTime: 2020-11-04 20:25:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:36
	% EndTime: 2020-11-04 20:25:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t39 = cos(qJ(1));
	t38 = sin(qJ(1));
	t1 = [t39, -t38, 0, 0; t38, t39, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:36
	% EndTime: 2020-11-04 20:25:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t41 = cos(qJ(1));
	t40 = sin(qJ(1));
	t1 = [0, -t41, t40, t41 * pkin(1) + t40 * qJ(2) + 0; 0, -t40, -t41, t40 * pkin(1) - t41 * qJ(2) + 0; 1, 0, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:36
	% EndTime: 2020-11-04 20:25:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t46 = pkin(1) + pkin(6);
	t45 = cos(qJ(1));
	t44 = cos(qJ(3));
	t43 = sin(qJ(1));
	t42 = sin(qJ(3));
	t1 = [t43 * t42, t43 * t44, t45, t43 * qJ(2) + t46 * t45 + 0; -t45 * t42, -t45 * t44, t43, -t45 * qJ(2) + t46 * t43 + 0; t44, -t42, 0, pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:36
	% EndTime: 2020-11-04 20:25:36
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (20->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->13)
	t48 = sin(qJ(4));
	t50 = sin(qJ(1));
	t58 = t50 * t48;
	t51 = cos(qJ(4));
	t57 = t50 * t51;
	t53 = cos(qJ(1));
	t56 = t53 * t48;
	t55 = t53 * t51;
	t54 = pkin(1) + pkin(6);
	t52 = cos(qJ(3));
	t49 = sin(qJ(3));
	t47 = -t49 * pkin(3) + t52 * pkin(7) - qJ(2);
	t1 = [t49 * t57 + t56, -t49 * t58 + t55, -t50 * t52, -t47 * t50 + t54 * t53 + 0; -t49 * t55 + t58, t49 * t56 + t57, t53 * t52, t47 * t53 + t54 * t50 + 0; t52 * t51, -t52 * t48, t49, t52 * pkin(3) + t49 * pkin(7) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:25:36
	% EndTime: 2020-11-04 20:25:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (28->20), mult. (31->22), div. (0->0), fcn. (44->6), ass. (0->15)
	t62 = sin(qJ(4));
	t64 = sin(qJ(1));
	t72 = t64 * t62;
	t65 = cos(qJ(4));
	t71 = t64 * t65;
	t67 = cos(qJ(1));
	t70 = t67 * t62;
	t69 = t67 * t65;
	t60 = t65 * pkin(4) + pkin(3);
	t61 = -qJ(5) - pkin(7);
	t63 = sin(qJ(3));
	t66 = cos(qJ(3));
	t68 = t60 * t63 + t61 * t66 + qJ(2);
	t59 = t62 * pkin(4) + pkin(1) + pkin(6);
	t1 = [t63 * t71 + t70, -t63 * t72 + t69, -t64 * t66, t59 * t67 + t68 * t64 + 0; -t63 * t69 + t72, t63 * t70 + t71, t67 * t66, t59 * t64 - t68 * t67 + 0; t66 * t65, -t66 * t62, t63, t66 * t60 - t61 * t63 + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end