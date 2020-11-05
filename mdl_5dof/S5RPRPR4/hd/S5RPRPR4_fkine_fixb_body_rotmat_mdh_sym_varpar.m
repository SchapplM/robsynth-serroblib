% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR4 (for one body)
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

function Tc_mdh = S5RPRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:19:04
	% EndTime: 2020-11-04 20:19:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:19:04
	% EndTime: 2020-11-04 20:19:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t37 = cos(qJ(1));
	t36 = sin(qJ(1));
	t1 = [0, 0, 1, pkin(5) + 0; t36, t37, 0, 0; -t37, t36, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:19:04
	% EndTime: 2020-11-04 20:19:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t40 = qJ(1) + pkin(8);
	t39 = cos(t40);
	t38 = sin(t40);
	t1 = [0, 0, 1, qJ(2) + pkin(5) + 0; t38, t39, 0, sin(qJ(1)) * pkin(1) + 0; -t39, t38, 0, -cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:19:04
	% EndTime: 2020-11-04 20:19:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->13), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t45 = cos(qJ(3));
	t44 = sin(qJ(3));
	t43 = qJ(1) + pkin(8);
	t42 = cos(t43);
	t41 = sin(t43);
	t1 = [t44, t45, 0, qJ(2) + pkin(5) + 0; t41 * t45, -t41 * t44, -t42, t41 * pkin(2) - t42 * pkin(6) + sin(qJ(1)) * pkin(1) + 0; -t42 * t45, t42 * t44, -t41, -t42 * pkin(2) - t41 * pkin(6) - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:19:04
	% EndTime: 2020-11-04 20:19:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->17), mult. (13->12), div. (0->0), fcn. (21->8), ass. (0->9)
	t53 = -qJ(4) - pkin(6);
	t52 = qJ(1) + pkin(8);
	t51 = qJ(3) + pkin(9);
	t50 = cos(t52);
	t49 = cos(t51);
	t48 = sin(t52);
	t47 = sin(t51);
	t46 = cos(qJ(3)) * pkin(3) + pkin(2);
	t1 = [t47, t49, 0, sin(qJ(3)) * pkin(3) + qJ(2) + pkin(5) + 0; t48 * t49, -t48 * t47, -t50, t48 * t46 + t50 * t53 + sin(qJ(1)) * pkin(1) + 0; -t50 * t49, t50 * t47, -t48, -t50 * t46 + t48 * t53 - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:19:04
	% EndTime: 2020-11-04 20:19:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (47->20), mult. (16->14), div. (0->0), fcn. (24->10), ass. (0->10)
	t61 = qJ(3) + pkin(9);
	t62 = qJ(1) + pkin(8);
	t60 = -pkin(7) - qJ(4) - pkin(6);
	t59 = qJ(5) + t61;
	t58 = cos(t62);
	t57 = sin(t62);
	t56 = cos(t59);
	t55 = sin(t59);
	t54 = -cos(qJ(3)) * pkin(3) - pkin(4) * cos(t61) - pkin(2);
	t1 = [t55, t56, 0, pkin(4) * sin(t61) + sin(qJ(3)) * pkin(3) + qJ(2) + pkin(5) + 0; t57 * t56, -t57 * t55, -t58, -t57 * t54 + t58 * t60 + sin(qJ(1)) * pkin(1) + 0; -t58 * t56, t58 * t55, -t57, t54 * t58 + t57 * t60 - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end