% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:26
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:26:42
	% EndTime: 2020-11-04 20:26:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:26:42
	% EndTime: 2020-11-04 20:26:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t36 = cos(qJ(1));
	t35 = sin(qJ(1));
	t1 = [0, 0, 1, pkin(5) + 0; t35, t36, 0, 0; -t36, t35, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:26:42
	% EndTime: 2020-11-04 20:26:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t39 = qJ(1) + pkin(9);
	t38 = cos(t39);
	t37 = sin(t39);
	t1 = [0, 0, 1, qJ(2) + pkin(5) + 0; t37, t38, 0, sin(qJ(1)) * pkin(1) + 0; -t38, t37, 0, -cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:26:42
	% EndTime: 2020-11-04 20:26:42
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (18->10), mult. (4->4), div. (0->0), fcn. (8->6), ass. (0->5)
	t43 = qJ(1) + pkin(9);
	t42 = qJ(3) + t43;
	t41 = cos(t42);
	t40 = sin(t42);
	t1 = [0, 0, 1, pkin(6) + qJ(2) + pkin(5) + 0; t40, t41, 0, pkin(2) * sin(t43) + sin(qJ(1)) * pkin(1) + 0; -t41, t40, 0, -pkin(2) * cos(t43) - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:26:42
	% EndTime: 2020-11-04 20:26:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (37->17), mult. (12->12), div. (0->0), fcn. (20->8), ass. (0->7)
	t47 = qJ(1) + pkin(9);
	t49 = cos(qJ(4));
	t48 = sin(qJ(4));
	t46 = qJ(3) + t47;
	t45 = cos(t46);
	t44 = sin(t46);
	t1 = [t48, t49, 0, pkin(6) + qJ(2) + pkin(5) + 0; t44 * t49, -t44 * t48, -t45, t44 * pkin(3) - t45 * pkin(7) + pkin(2) * sin(t47) + sin(qJ(1)) * pkin(1) + 0; -t45 * t49, t45 * t48, -t44, -t45 * pkin(3) - t44 * pkin(7) - pkin(2) * cos(t47) - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:26:42
	% EndTime: 2020-11-04 20:26:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (48->21), mult. (15->14), div. (0->0), fcn. (23->10), ass. (0->10)
	t56 = qJ(1) + pkin(9);
	t58 = -pkin(8) - pkin(7);
	t57 = qJ(4) + qJ(5);
	t55 = cos(t57);
	t54 = sin(t57);
	t53 = qJ(3) + t56;
	t52 = cos(qJ(4)) * pkin(4) + pkin(3);
	t51 = cos(t53);
	t50 = sin(t53);
	t1 = [t54, t55, 0, sin(qJ(4)) * pkin(4) + pkin(6) + qJ(2) + pkin(5) + 0; t50 * t55, -t50 * t54, -t51, t50 * t52 + t51 * t58 + pkin(2) * sin(t56) + sin(qJ(1)) * pkin(1) + 0; -t51 * t55, t51 * t54, -t50, -t51 * t52 + t50 * t58 - pkin(2) * cos(t56) - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end