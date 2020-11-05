% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRP3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:12
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:38
	% EndTime: 2020-11-04 20:12:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:38
	% EndTime: 2020-11-04 20:12:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t34 = cos(qJ(1));
	t33 = sin(qJ(1));
	t1 = [t34, -t33, 0, 0; t33, t34, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:38
	% EndTime: 2020-11-04 20:12:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t37 = qJ(1) + pkin(7);
	t36 = cos(t37);
	t35 = sin(t37);
	t1 = [t36, -t35, 0, cos(qJ(1)) * pkin(1) + 0; t35, t36, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:38
	% EndTime: 2020-11-04 20:12:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (6->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t40 = qJ(1) + pkin(7);
	t39 = cos(t40);
	t38 = sin(t40);
	t1 = [0, -t39, t38, t39 * pkin(2) + t38 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; 0, -t38, -t39, t38 * pkin(2) - t39 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; 1, 0, 0, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:38
	% EndTime: 2020-11-04 20:12:38
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (24->14), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->7)
	t46 = pkin(2) + pkin(6);
	t45 = cos(qJ(4));
	t44 = sin(qJ(4));
	t43 = qJ(1) + pkin(7);
	t42 = cos(t43);
	t41 = sin(t43);
	t1 = [t41 * t44, t41 * t45, t42, t46 * t42 + t41 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; -t42 * t44, -t42 * t45, t41, t46 * t41 - t42 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; t45, -t44, 0, pkin(3) + qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:12:38
	% EndTime: 2020-11-04 20:12:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->16), mult. (17->12), div. (0->0), fcn. (25->6), ass. (0->8)
	t54 = pkin(2) + qJ(5) + pkin(6);
	t51 = sin(qJ(4));
	t53 = pkin(4) * t51 + qJ(3);
	t52 = cos(qJ(4));
	t49 = qJ(1) + pkin(7);
	t48 = cos(t49);
	t47 = sin(t49);
	t1 = [t47 * t51, t47 * t52, t48, cos(qJ(1)) * pkin(1) + 0 + t54 * t48 + t53 * t47; -t48 * t51, -t48 * t52, t47, sin(qJ(1)) * pkin(1) + 0 - t53 * t48 + t54 * t47; t52, -t51, 0, t52 * pkin(4) + pkin(3) + pkin(5) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end