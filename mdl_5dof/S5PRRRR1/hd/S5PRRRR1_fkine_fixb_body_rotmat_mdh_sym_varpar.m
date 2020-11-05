% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:07
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [2x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:33
	% EndTime: 2020-11-04 20:07:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:33
	% EndTime: 2020-11-04 20:07:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:33
	% EndTime: 2020-11-04 20:07:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t32 = cos(qJ(2));
	t31 = sin(qJ(2));
	t1 = [t32, -t31, 0, pkin(1) + 0; 0, 0, -1, 0; t31, t32, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:33
	% EndTime: 2020-11-04 20:07:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (7->7), mult. (4->4), div. (0->0), fcn. (12->4), ass. (0->5)
	t36 = cos(qJ(2));
	t35 = cos(qJ(3));
	t34 = sin(qJ(2));
	t33 = sin(qJ(3));
	t1 = [t36 * t35, -t36 * t33, t34, pkin(1) + 0; -t33, -t35, 0, 0; t34 * t35, -t34 * t33, -t36, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:33
	% EndTime: 2020-11-04 20:07:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (16->11), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->7)
	t43 = pkin(2) * cos(qJ(3));
	t42 = cos(qJ(2));
	t40 = sin(qJ(2));
	t39 = qJ(3) + qJ(4);
	t38 = cos(t39);
	t37 = sin(t39);
	t1 = [t42 * t38, -t42 * t37, t40, t42 * t43 + pkin(1) + 0; -t37, -t38, 0, -sin(qJ(3)) * pkin(2) + 0; t40 * t38, -t40 * t37, -t42, t40 * t43 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:33
	% EndTime: 2020-11-04 20:07:33
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (19->11), mult. (21->16), div. (0->0), fcn. (34->8), ass. (0->13)
	t56 = pkin(2) * cos(qJ(3));
	t47 = sin(qJ(5));
	t48 = sin(qJ(2));
	t55 = t48 * t47;
	t49 = cos(qJ(5));
	t54 = t48 * t49;
	t51 = cos(qJ(2));
	t53 = t51 * t47;
	t52 = t51 * t49;
	t46 = qJ(3) + qJ(4);
	t45 = cos(t46);
	t44 = sin(t46);
	t1 = [t45 * t52 + t55, -t45 * t53 + t54, t51 * t44, t51 * t56 + pkin(1) + 0; -t44 * t49, t44 * t47, t45, -sin(qJ(3)) * pkin(2) + 0; t45 * t54 - t53, -t45 * t55 - t52, t48 * t44, t48 * t56 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end