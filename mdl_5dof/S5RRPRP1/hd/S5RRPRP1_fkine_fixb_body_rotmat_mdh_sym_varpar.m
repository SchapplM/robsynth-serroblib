% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRP1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:32
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPRP1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:32:30
	% EndTime: 2020-11-04 20:32:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:32:30
	% EndTime: 2020-11-04 20:32:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t35 = cos(qJ(1));
	t34 = sin(qJ(1));
	t1 = [0, 0, 1, pkin(5) + 0; t34, t35, 0, 0; -t35, t34, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:32:30
	% EndTime: 2020-11-04 20:32:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t38 = qJ(1) + qJ(2);
	t37 = cos(t38);
	t36 = sin(t38);
	t1 = [0, 0, 1, pkin(6) + pkin(5) + 0; t36, t37, 0, sin(qJ(1)) * pkin(1) + 0; -t37, t36, 0, -cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:32:30
	% EndTime: 2020-11-04 20:32:30
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (18->10), mult. (4->4), div. (0->0), fcn. (8->6), ass. (0->5)
	t42 = qJ(1) + qJ(2);
	t41 = pkin(8) + t42;
	t40 = cos(t41);
	t39 = sin(t41);
	t1 = [0, 0, 1, qJ(3) + pkin(6) + pkin(5) + 0; t39, t40, 0, pkin(2) * sin(t42) + sin(qJ(1)) * pkin(1) + 0; -t40, t39, 0, -pkin(2) * cos(t42) - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:32:30
	% EndTime: 2020-11-04 20:32:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (37->17), mult. (12->12), div. (0->0), fcn. (20->8), ass. (0->7)
	t46 = qJ(1) + qJ(2);
	t48 = cos(qJ(4));
	t47 = sin(qJ(4));
	t45 = pkin(8) + t46;
	t44 = cos(t45);
	t43 = sin(t45);
	t1 = [t47, t48, 0, qJ(3) + pkin(6) + pkin(5) + 0; t43 * t48, -t43 * t47, -t44, t43 * pkin(3) - t44 * pkin(7) + pkin(2) * sin(t46) + sin(qJ(1)) * pkin(1) + 0; -t44 * t48, t44 * t47, -t43, -t44 * pkin(3) - t43 * pkin(7) - pkin(2) * cos(t46) - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:32:30
	% EndTime: 2020-11-04 20:32:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (42->20), mult. (15->14), div. (0->0), fcn. (23->8), ass. (0->9)
	t53 = qJ(1) + qJ(2);
	t56 = cos(qJ(4));
	t55 = sin(qJ(4));
	t54 = -qJ(5) - pkin(7);
	t52 = pkin(8) + t53;
	t51 = t56 * pkin(4) + pkin(3);
	t50 = cos(t52);
	t49 = sin(t52);
	t1 = [t55, t56, 0, t55 * pkin(4) + pkin(5) + pkin(6) + qJ(3) + 0; t49 * t56, -t49 * t55, -t50, t49 * t51 + t50 * t54 + pkin(2) * sin(t53) + sin(qJ(1)) * pkin(1) + 0; -t50 * t56, t50 * t55, -t49, -t50 * t51 + t49 * t54 - pkin(2) * cos(t53) - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end