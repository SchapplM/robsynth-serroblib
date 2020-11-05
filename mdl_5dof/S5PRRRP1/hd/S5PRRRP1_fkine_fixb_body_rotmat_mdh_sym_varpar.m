% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRP1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:05
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:05:20
	% EndTime: 2020-11-04 20:05:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:05:20
	% EndTime: 2020-11-04 20:05:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t37 = cos(pkin(8));
	t36 = sin(pkin(8));
	t1 = [t37, -t36, 0, 0; t36, t37, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:05:20
	% EndTime: 2020-11-04 20:05:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t40 = pkin(8) + qJ(2);
	t39 = cos(t40);
	t38 = sin(t40);
	t1 = [t39, -t38, 0, cos(pkin(8)) * pkin(1) + 0; t38, t39, 0, sin(pkin(8)) * pkin(1) + 0; 0, 0, 1, pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:05:20
	% EndTime: 2020-11-04 20:05:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->10), mult. (4->4), div. (0->0), fcn. (8->6), ass. (0->5)
	t44 = pkin(8) + qJ(2);
	t43 = qJ(3) + t44;
	t42 = cos(t43);
	t41 = sin(t43);
	t1 = [t42, -t41, 0, pkin(2) * cos(t44) + cos(pkin(8)) * pkin(1) + 0; t41, t42, 0, pkin(2) * sin(t44) + sin(pkin(8)) * pkin(1) + 0; 0, 0, 1, pkin(6) + pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:05:20
	% EndTime: 2020-11-04 20:05:20
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (36->16), mult. (12->12), div. (0->0), fcn. (20->8), ass. (0->7)
	t48 = pkin(8) + qJ(2);
	t50 = cos(qJ(4));
	t49 = sin(qJ(4));
	t47 = qJ(3) + t48;
	t46 = cos(t47);
	t45 = sin(t47);
	t1 = [t46 * t50, -t46 * t49, t45, t46 * pkin(3) + t45 * pkin(7) + pkin(2) * cos(t48) + cos(pkin(8)) * pkin(1) + 0; t45 * t50, -t45 * t49, -t46, t45 * pkin(3) - t46 * pkin(7) + pkin(2) * sin(t48) + sin(pkin(8)) * pkin(1) + 0; t49, t50, 0, pkin(6) + pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:05:20
	% EndTime: 2020-11-04 20:05:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (41->19), mult. (15->14), div. (0->0), fcn. (23->8), ass. (0->9)
	t55 = pkin(8) + qJ(2);
	t58 = cos(qJ(4));
	t57 = sin(qJ(4));
	t56 = qJ(5) + pkin(7);
	t54 = qJ(3) + t55;
	t53 = t58 * pkin(4) + pkin(3);
	t52 = cos(t54);
	t51 = sin(t54);
	t1 = [t52 * t58, -t52 * t57, t51, t52 * t53 + t56 * t51 + cos(pkin(8)) * pkin(1) + pkin(2) * cos(t55) + 0; t51 * t58, -t51 * t57, -t52, t51 * t53 - t52 * t56 + pkin(2) * sin(t55) + sin(pkin(8)) * pkin(1) + 0; t57, t58, 0, t57 * pkin(4) + pkin(5) + pkin(6) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end