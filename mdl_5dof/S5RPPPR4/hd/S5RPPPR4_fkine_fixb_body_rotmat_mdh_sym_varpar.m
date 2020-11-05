% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPPR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:11
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:16
	% EndTime: 2020-11-04 20:11:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:16
	% EndTime: 2020-11-04 20:11:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t36 = cos(qJ(1));
	t35 = sin(qJ(1));
	t1 = [t36, -t35, 0, 0; t35, t36, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:16
	% EndTime: 2020-11-04 20:11:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t39 = qJ(1) + pkin(7);
	t38 = cos(t39);
	t37 = sin(t39);
	t1 = [t38, -t37, 0, cos(qJ(1)) * pkin(1) + 0; t37, t38, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:16
	% EndTime: 2020-11-04 20:11:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (6->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t42 = qJ(1) + pkin(7);
	t41 = cos(t42);
	t40 = sin(t42);
	t1 = [0, -t41, t40, t41 * pkin(2) + t40 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; 0, -t40, -t41, t40 * pkin(2) - t41 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; 1, 0, 0, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:16
	% EndTime: 2020-11-04 20:11:16
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (24->14), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->7)
	t48 = pkin(2) + qJ(4);
	t47 = cos(pkin(8));
	t46 = sin(pkin(8));
	t45 = qJ(1) + pkin(7);
	t44 = cos(t45);
	t43 = sin(t45);
	t1 = [t43 * t46, t43 * t47, t44, t48 * t44 + t43 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; -t44 * t46, -t44 * t47, t43, t48 * t43 - t44 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; t47, -t46, 0, pkin(3) + qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:11:16
	% EndTime: 2020-11-04 20:11:16
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (39->17), mult. (17->12), div. (0->0), fcn. (25->8), ass. (0->9)
	t58 = pkin(2) + pkin(6) + qJ(4);
	t57 = sin(pkin(8)) * pkin(4) + qJ(3);
	t54 = qJ(1) + pkin(7);
	t53 = pkin(8) + qJ(5);
	t52 = cos(t54);
	t51 = cos(t53);
	t50 = sin(t54);
	t49 = sin(t53);
	t1 = [t50 * t49, t50 * t51, t52, cos(qJ(1)) * pkin(1) + 0 + t58 * t52 + t57 * t50; -t52 * t49, -t52 * t51, t50, sin(qJ(1)) * pkin(1) + 0 - t57 * t52 + t58 * t50; t51, -t49, 0, cos(pkin(8)) * pkin(4) + pkin(3) + qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end