% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR11 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:16
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPRR11_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR11_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:24
	% EndTime: 2020-11-04 20:16:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:24
	% EndTime: 2020-11-04 20:16:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t35 = cos(qJ(1));
	t34 = sin(qJ(1));
	t1 = [t35, -t34, 0, 0; t34, t35, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:24
	% EndTime: 2020-11-04 20:16:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t37 = cos(qJ(1));
	t36 = sin(qJ(1));
	t1 = [0, -t37, t36, t37 * pkin(1) + t36 * qJ(2) + 0; 0, -t36, -t37, t36 * pkin(1) - t37 * qJ(2) + 0; 1, 0, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:24
	% EndTime: 2020-11-04 20:16:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->4)
	t40 = cos(qJ(1));
	t39 = sin(qJ(1));
	t38 = pkin(1) + qJ(3);
	t1 = [0, t39, t40, t39 * qJ(2) + t38 * t40 + 0; 0, -t40, t39, -t40 * qJ(2) + t38 * t39 + 0; 1, 0, 0, pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:24
	% EndTime: 2020-11-04 20:16:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->11), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->7)
	t46 = cos(qJ(1));
	t45 = cos(qJ(4));
	t44 = sin(qJ(1));
	t43 = sin(qJ(4));
	t42 = pkin(1) + qJ(3);
	t41 = pkin(6) - qJ(2);
	t1 = [t46 * t43, t46 * t45, -t44, -t41 * t44 + t42 * t46 + 0; t44 * t43, t44 * t45, t46, t41 * t46 + t42 * t44 + 0; t45, -t43, 0, pkin(3) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:16:24
	% EndTime: 2020-11-04 20:16:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (24->20), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->13)
	t49 = sin(qJ(5));
	t51 = sin(qJ(1));
	t58 = t51 * t49;
	t52 = cos(qJ(5));
	t57 = t51 * t52;
	t54 = cos(qJ(1));
	t56 = t54 * t49;
	t55 = t54 * t52;
	t53 = cos(qJ(4));
	t50 = sin(qJ(4));
	t48 = pkin(6) - qJ(2);
	t47 = t50 * pkin(4) - t53 * pkin(7) + pkin(1) + qJ(3);
	t1 = [t50 * t55 - t58, -t50 * t56 - t57, -t54 * t53, t47 * t54 - t48 * t51 + 0; t50 * t57 + t56, -t50 * t58 + t55, -t51 * t53, t47 * t51 + t48 * t54 + 0; t53 * t52, -t53 * t49, t50, t53 * pkin(4) + t50 * pkin(7) + pkin(2) + pkin(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end