% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:59
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:59:25
	% EndTime: 2020-11-04 19:59:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:59:25
	% EndTime: 2020-11-04 19:59:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t38 = cos(pkin(8));
	t37 = sin(pkin(8));
	t1 = [t38, -t37, 0, 0; t37, t38, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:59:25
	% EndTime: 2020-11-04 19:59:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t41 = pkin(8) + qJ(2);
	t40 = cos(t41);
	t39 = sin(t41);
	t1 = [t40, -t39, 0, cos(pkin(8)) * pkin(1) + 0; t39, t40, 0, sin(pkin(8)) * pkin(1) + 0; 0, 0, 1, pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:59:25
	% EndTime: 2020-11-04 19:59:25
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t46 = cos(pkin(9));
	t45 = sin(pkin(9));
	t44 = pkin(8) + qJ(2);
	t43 = cos(t44);
	t42 = sin(t44);
	t1 = [t43 * t46, -t43 * t45, t42, t43 * pkin(2) + t42 * qJ(3) + cos(pkin(8)) * pkin(1) + 0; t42 * t46, -t42 * t45, -t43, t42 * pkin(2) - t43 * qJ(3) + sin(pkin(8)) * pkin(1) + 0; t45, t46, 0, pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:59:25
	% EndTime: 2020-11-04 19:59:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (32->16), mult. (13->12), div. (0->0), fcn. (21->8), ass. (0->9)
	t54 = -pkin(6) - qJ(3);
	t53 = pkin(8) + qJ(2);
	t52 = pkin(9) + qJ(4);
	t51 = cos(t53);
	t50 = cos(t52);
	t49 = sin(t53);
	t48 = sin(t52);
	t47 = cos(pkin(9)) * pkin(3) + pkin(2);
	t1 = [t51 * t50, -t51 * t48, t49, t51 * t47 - t49 * t54 + cos(pkin(8)) * pkin(1) + 0; t49 * t50, -t49 * t48, -t51, t49 * t47 + t51 * t54 + sin(pkin(8)) * pkin(1) + 0; t48, t50, 0, sin(pkin(9)) * pkin(3) + pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:59:25
	% EndTime: 2020-11-04 19:59:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (46->19), mult. (16->14), div. (0->0), fcn. (24->10), ass. (0->10)
	t62 = pkin(9) + qJ(4);
	t63 = pkin(8) + qJ(2);
	t61 = qJ(3) + pkin(6) + pkin(7);
	t60 = qJ(5) + t62;
	t59 = cos(t63);
	t58 = sin(t63);
	t57 = cos(t60);
	t56 = sin(t60);
	t55 = pkin(4) * cos(t62) + cos(pkin(9)) * pkin(3) + pkin(2);
	t1 = [t59 * t57, -t59 * t56, t58, t59 * t55 + t61 * t58 + cos(pkin(8)) * pkin(1) + 0; t58 * t57, -t58 * t56, -t59, t58 * t55 - t59 * t61 + sin(pkin(8)) * pkin(1) + 0; t56, t57, 0, pkin(4) * sin(t62) + sin(pkin(9)) * pkin(3) + pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end