% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP3 (for one body)
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

function Tc_mdh = S5RPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:09
	% EndTime: 2020-11-04 20:23:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:09
	% EndTime: 2020-11-04 20:23:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t36 = cos(qJ(1));
	t35 = sin(qJ(1));
	t1 = [0, 0, 1, pkin(5) + 0; t35, t36, 0, 0; -t36, t35, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:09
	% EndTime: 2020-11-04 20:23:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t39 = qJ(1) + pkin(8);
	t38 = cos(t39);
	t37 = sin(t39);
	t1 = [0, 0, 1, qJ(2) + pkin(5) + 0; t37, t38, 0, sin(qJ(1)) * pkin(1) + 0; -t38, t37, 0, -cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:09
	% EndTime: 2020-11-04 20:23:09
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (22->13), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t44 = cos(qJ(3));
	t43 = sin(qJ(3));
	t42 = qJ(1) + pkin(8);
	t41 = cos(t42);
	t40 = sin(t42);
	t1 = [t43, t44, 0, qJ(2) + pkin(5) + 0; t40 * t44, -t40 * t43, -t41, t40 * pkin(2) - t41 * pkin(6) + sin(qJ(1)) * pkin(1) + 0; -t41 * t44, t41 * t43, -t40, -t41 * pkin(2) - t40 * pkin(6) - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:09
	% EndTime: 2020-11-04 20:23:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->17), mult. (13->12), div. (0->0), fcn. (21->8), ass. (0->9)
	t52 = -pkin(7) - pkin(6);
	t51 = qJ(3) + qJ(4);
	t50 = qJ(1) + pkin(8);
	t49 = cos(t51);
	t48 = sin(t51);
	t47 = cos(t50);
	t46 = sin(t50);
	t45 = cos(qJ(3)) * pkin(3) + pkin(2);
	t1 = [t48, t49, 0, sin(qJ(3)) * pkin(3) + qJ(2) + pkin(5) + 0; t46 * t49, -t46 * t48, -t47, t46 * t45 + t47 * t52 + sin(qJ(1)) * pkin(1) + 0; -t47 * t49, t47 * t48, -t46, -t47 * t45 + t46 * t52 - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:09
	% EndTime: 2020-11-04 20:23:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (41->19), mult. (16->14), div. (0->0), fcn. (24->8), ass. (0->9)
	t60 = qJ(3) + qJ(4);
	t59 = qJ(1) + pkin(8);
	t58 = qJ(5) + pkin(7) + pkin(6);
	t57 = cos(t60);
	t56 = sin(t60);
	t55 = cos(t59);
	t54 = sin(t59);
	t53 = pkin(4) * t57 + cos(qJ(3)) * pkin(3) + pkin(2);
	t1 = [t56, t57, 0, pkin(4) * t56 + sin(qJ(3)) * pkin(3) + qJ(2) + pkin(5) + 0; t54 * t57, -t54 * t56, -t55, t54 * t53 - t55 * t58 + sin(qJ(1)) * pkin(1) + 0; -t55 * t57, t55 * t56, -t54, -t55 * t53 - t54 * t58 - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end