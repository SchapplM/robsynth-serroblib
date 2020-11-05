% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR16 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:22
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRPR16_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR16_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:21
	% EndTime: 2020-11-04 20:22:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:21
	% EndTime: 2020-11-04 20:22:21
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t33 = cos(qJ(1));
	t32 = sin(qJ(1));
	t1 = [t33, -t32, 0, 0; t32, t33, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:21
	% EndTime: 2020-11-04 20:22:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t35 = cos(qJ(1));
	t34 = sin(qJ(1));
	t1 = [0, -t35, t34, t35 * pkin(1) + t34 * qJ(2) + 0; 0, -t34, -t35, t34 * pkin(1) - t35 * qJ(2) + 0; 1, 0, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:21
	% EndTime: 2020-11-04 20:22:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t40 = pkin(1) + pkin(6);
	t39 = cos(qJ(1));
	t38 = cos(qJ(3));
	t37 = sin(qJ(1));
	t36 = sin(qJ(3));
	t1 = [t37 * t36, t37 * t38, t39, t37 * qJ(2) + t40 * t39 + 0; -t39 * t36, -t39 * t38, t37, -t39 * qJ(2) + t40 * t37 + 0; t38, -t36, 0, pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:21
	% EndTime: 2020-11-04 20:22:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (17->14), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->7)
	t46 = pkin(1) + pkin(6);
	t45 = cos(qJ(1));
	t44 = cos(qJ(3));
	t43 = sin(qJ(1));
	t42 = sin(qJ(3));
	t41 = -t42 * pkin(3) + t44 * qJ(4) - qJ(2);
	t1 = [t45, -t43 * t42, -t43 * t44, -t41 * t43 + t46 * t45 + 0; t43, t45 * t42, t45 * t44, t41 * t45 + t46 * t43 + 0; 0, -t44, t42, t44 * pkin(3) + t42 * qJ(4) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:21
	% EndTime: 2020-11-04 20:22:21
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (24->17), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->13)
	t50 = sin(qJ(1));
	t52 = cos(qJ(3));
	t58 = t50 * t52;
	t48 = sin(qJ(5));
	t53 = cos(qJ(1));
	t57 = t53 * t48;
	t51 = cos(qJ(5));
	t56 = t53 * t51;
	t49 = sin(qJ(3));
	t54 = pkin(3) + pkin(7);
	t55 = t52 * qJ(4) - t54 * t49 - qJ(2);
	t47 = pkin(1) + pkin(4) + pkin(6);
	t1 = [-t48 * t58 + t56, -t51 * t58 - t57, t50 * t49, t47 * t53 - t55 * t50 + 0; t50 * t51 + t52 * t57, -t50 * t48 + t52 * t56, -t53 * t49, t47 * t50 + t55 * t53 + 0; t49 * t48, t49 * t51, t52, t49 * qJ(4) + t54 * t52 + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end