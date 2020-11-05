% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:36
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:36:00
	% EndTime: 2020-11-04 20:36:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:36:00
	% EndTime: 2020-11-04 20:36:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t38 = cos(qJ(1));
	t37 = sin(qJ(1));
	t1 = [0, 0, 1, pkin(5) + 0; t37, t38, 0, 0; -t38, t37, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:36:00
	% EndTime: 2020-11-04 20:36:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t41 = qJ(1) + qJ(2);
	t40 = cos(t41);
	t39 = sin(t41);
	t1 = [0, 0, 1, pkin(6) + pkin(5) + 0; t39, t40, 0, sin(qJ(1)) * pkin(1) + 0; -t40, t39, 0, -cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:36:00
	% EndTime: 2020-11-04 20:36:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->10), mult. (4->4), div. (0->0), fcn. (8->6), ass. (0->5)
	t45 = qJ(1) + qJ(2);
	t44 = pkin(9) + t45;
	t43 = cos(t44);
	t42 = sin(t44);
	t1 = [0, 0, 1, qJ(3) + pkin(6) + pkin(5) + 0; t42, t43, 0, pkin(2) * sin(t45) + sin(qJ(1)) * pkin(1) + 0; -t43, t42, 0, -pkin(2) * cos(t45) - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:36:00
	% EndTime: 2020-11-04 20:36:00
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->14), mult. (6->6), div. (0->0), fcn. (10->8), ass. (0->6)
	t50 = qJ(1) + qJ(2);
	t49 = pkin(9) + t50;
	t48 = qJ(4) + t49;
	t47 = cos(t48);
	t46 = sin(t48);
	t1 = [0, 0, 1, pkin(7) + qJ(3) + pkin(6) + pkin(5) + 0; t46, t47, 0, pkin(3) * sin(t49) + pkin(2) * sin(t50) + sin(qJ(1)) * pkin(1) + 0; -t47, t46, 0, -pkin(3) * cos(t49) - pkin(2) * cos(t50) - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:36:00
	% EndTime: 2020-11-04 20:36:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (54->21), mult. (14->14), div. (0->0), fcn. (22->10), ass. (0->8)
	t55 = qJ(1) + qJ(2);
	t54 = pkin(9) + t55;
	t57 = cos(qJ(5));
	t56 = sin(qJ(5));
	t53 = qJ(4) + t54;
	t52 = cos(t53);
	t51 = sin(t53);
	t1 = [t56, t57, 0, pkin(7) + qJ(3) + pkin(6) + pkin(5) + 0; t51 * t57, -t51 * t56, -t52, t51 * pkin(4) - t52 * pkin(8) + pkin(3) * sin(t54) + pkin(2) * sin(t55) + sin(qJ(1)) * pkin(1) + 0; -t52 * t57, t52 * t56, -t51, -t52 * pkin(4) - t51 * pkin(8) - pkin(3) * cos(t54) - pkin(2) * cos(t55) - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end