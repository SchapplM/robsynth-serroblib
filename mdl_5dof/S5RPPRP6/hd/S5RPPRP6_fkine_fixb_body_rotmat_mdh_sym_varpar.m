% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRP6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:13
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:25
	% EndTime: 2020-11-04 20:13:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:25
	% EndTime: 2020-11-04 20:13:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t36 = cos(qJ(1));
	t35 = sin(qJ(1));
	t1 = [t36, -t35, 0, 0; t35, t36, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:25
	% EndTime: 2020-11-04 20:13:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t38 = cos(qJ(1));
	t37 = sin(qJ(1));
	t1 = [0, -t38, t37, t38 * pkin(1) + t37 * qJ(2) + 0; 0, -t37, -t38, t37 * pkin(1) - t38 * qJ(2) + 0; 1, 0, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:25
	% EndTime: 2020-11-04 20:13:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t43 = cos(qJ(1));
	t42 = sin(qJ(1));
	t41 = pkin(1) + qJ(3);
	t40 = cos(pkin(7));
	t39 = sin(pkin(7));
	t1 = [t42 * t39, t42 * t40, t43, t42 * qJ(2) + t41 * t43 + 0; -t43 * t39, -t43 * t40, t42, -t43 * qJ(2) + t41 * t42 + 0; t40, -t39, 0, pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:25
	% EndTime: 2020-11-04 20:13:25
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (22->13), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t50 = cos(qJ(1));
	t49 = sin(qJ(1));
	t48 = pkin(7) + qJ(4);
	t47 = pkin(1) + pkin(6) + qJ(3);
	t46 = cos(t48);
	t45 = sin(t48);
	t44 = sin(pkin(7)) * pkin(3) + qJ(2);
	t1 = [t49 * t45, t49 * t46, t50, t44 * t49 + t47 * t50 + 0; -t50 * t45, -t50 * t46, t49, -t44 * t50 + t47 * t49 + 0; t46, -t45, 0, cos(pkin(7)) * pkin(3) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:25
	% EndTime: 2020-11-04 20:13:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->16), mult. (21->14), div. (0->0), fcn. (29->6), ass. (0->8)
	t55 = pkin(7) + qJ(4);
	t52 = sin(t55);
	t53 = cos(t55);
	t58 = pkin(4) * t52 - qJ(5) * t53 + sin(pkin(7)) * pkin(3) + qJ(2);
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t54 = pkin(1) + pkin(6) + qJ(3);
	t1 = [t56 * t52, t57, -t56 * t53, t54 * t57 + t58 * t56 + 0; -t57 * t52, t56, t57 * t53, t54 * t56 - t58 * t57 + 0; t53, 0, t52, t53 * pkin(4) + t52 * qJ(5) + cos(pkin(7)) * pkin(3) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end