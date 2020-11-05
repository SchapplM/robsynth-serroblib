% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR15 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:22
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRPR15_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR15_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:05
	% EndTime: 2020-11-04 20:22:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:05
	% EndTime: 2020-11-04 20:22:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t36 = cos(qJ(1));
	t35 = sin(qJ(1));
	t1 = [t36, -t35, 0, 0; t35, t36, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:05
	% EndTime: 2020-11-04 20:22:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t38 = cos(qJ(1));
	t37 = sin(qJ(1));
	t1 = [0, -t38, t37, t38 * pkin(1) + t37 * qJ(2) + 0; 0, -t37, -t38, t37 * pkin(1) - t38 * qJ(2) + 0; 1, 0, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:05
	% EndTime: 2020-11-04 20:22:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t43 = pkin(1) + pkin(6);
	t42 = cos(qJ(1));
	t41 = cos(qJ(3));
	t40 = sin(qJ(1));
	t39 = sin(qJ(3));
	t1 = [t40 * t39, t40 * t41, t42, t40 * qJ(2) + t43 * t42 + 0; -t42 * t39, -t42 * t41, t40, -t42 * qJ(2) + t43 * t40 + 0; t41, -t39, 0, pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:05
	% EndTime: 2020-11-04 20:22:05
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (20->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->13)
	t45 = sin(pkin(8));
	t48 = sin(qJ(1));
	t55 = t48 * t45;
	t46 = cos(pkin(8));
	t54 = t48 * t46;
	t50 = cos(qJ(1));
	t53 = t50 * t45;
	t52 = t50 * t46;
	t51 = pkin(1) + pkin(6);
	t49 = cos(qJ(3));
	t47 = sin(qJ(3));
	t44 = -t47 * pkin(3) + t49 * qJ(4) - qJ(2);
	t1 = [t47 * t54 + t53, -t47 * t55 + t52, -t48 * t49, -t44 * t48 + t51 * t50 + 0; -t47 * t52 + t55, t47 * t53 + t54, t50 * t49, t44 * t50 + t51 * t48 + 0; t49 * t46, -t49 * t45, t47, t49 * pkin(3) + t47 * qJ(4) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:05
	% EndTime: 2020-11-04 20:22:05
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (38->21), mult. (31->22), div. (0->0), fcn. (44->8), ass. (0->16)
	t60 = pkin(8) + qJ(5);
	t58 = sin(t60);
	t63 = sin(qJ(1));
	t70 = t63 * t58;
	t59 = cos(t60);
	t69 = t63 * t59;
	t65 = cos(qJ(1));
	t68 = t65 * t58;
	t67 = t65 * t59;
	t57 = cos(pkin(8)) * pkin(4) + pkin(3);
	t61 = pkin(7) + qJ(4);
	t62 = sin(qJ(3));
	t64 = cos(qJ(3));
	t66 = t57 * t62 - t61 * t64 + qJ(2);
	t56 = sin(pkin(8)) * pkin(4) + pkin(1) + pkin(6);
	t1 = [t62 * t69 + t68, -t62 * t70 + t67, -t63 * t64, t56 * t65 + t66 * t63 + 0; -t62 * t67 + t70, t62 * t68 + t69, t65 * t64, t56 * t63 - t66 * t65 + 0; t64 * t59, -t64 * t58, t62, t64 * t57 + t61 * t62 + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end