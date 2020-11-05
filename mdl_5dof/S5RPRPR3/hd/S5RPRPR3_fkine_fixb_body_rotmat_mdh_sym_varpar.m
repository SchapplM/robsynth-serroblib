% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:18
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:18:48
	% EndTime: 2020-11-04 20:18:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:18:48
	% EndTime: 2020-11-04 20:18:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t41 = cos(qJ(1));
	t40 = sin(qJ(1));
	t1 = [0, 0, 1, pkin(5) + 0; t40, t41, 0, 0; -t41, t40, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:18:48
	% EndTime: 2020-11-04 20:18:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t44 = qJ(1) + pkin(8);
	t43 = cos(t44);
	t42 = sin(t44);
	t1 = [0, 0, 1, qJ(2) + pkin(5) + 0; t42, t43, 0, sin(qJ(1)) * pkin(1) + 0; -t43, t42, 0, -cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:18:48
	% EndTime: 2020-11-04 20:18:48
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (18->10), mult. (4->4), div. (0->0), fcn. (8->6), ass. (0->5)
	t48 = qJ(1) + pkin(8);
	t47 = qJ(3) + t48;
	t46 = cos(t47);
	t45 = sin(t47);
	t1 = [0, 0, 1, pkin(6) + qJ(2) + pkin(5) + 0; t45, t46, 0, pkin(2) * sin(t48) + sin(qJ(1)) * pkin(1) + 0; -t46, t45, 0, -pkin(2) * cos(t48) - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:18:48
	% EndTime: 2020-11-04 20:18:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (37->17), mult. (12->12), div. (0->0), fcn. (20->8), ass. (0->7)
	t52 = qJ(1) + pkin(8);
	t54 = cos(pkin(9));
	t53 = sin(pkin(9));
	t51 = qJ(3) + t52;
	t50 = cos(t51);
	t49 = sin(t51);
	t1 = [t53, t54, 0, pkin(6) + qJ(2) + pkin(5) + 0; t49 * t54, -t49 * t53, -t50, t49 * pkin(3) - t50 * qJ(4) + pkin(2) * sin(t52) + sin(qJ(1)) * pkin(1) + 0; -t50 * t54, t50 * t53, -t49, -t50 * pkin(3) - t49 * qJ(4) - pkin(2) * cos(t52) - cos(qJ(1)) * pkin(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:18:48
	% EndTime: 2020-11-04 20:18:48
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (62->24), mult. (34->26), div. (0->0), fcn. (47->10), ass. (0->12)
	t60 = cos(pkin(9));
	t61 = sin(qJ(5));
	t65 = t60 * t61;
	t62 = cos(qJ(5));
	t64 = t60 * t62;
	t58 = qJ(1) + pkin(8);
	t59 = sin(pkin(9));
	t63 = pkin(4) * t60 + pkin(7) * t59 + pkin(3);
	t57 = qJ(3) + t58;
	t56 = cos(t57);
	t55 = sin(t57);
	t1 = [t59 * t62, -t59 * t61, -t60, t59 * pkin(4) - t60 * pkin(7) + pkin(5) + pkin(6) + qJ(2) + 0; t55 * t64 - t56 * t61, -t55 * t65 - t56 * t62, t55 * t59, pkin(2) * sin(t58) + sin(qJ(1)) * pkin(1) - t56 * qJ(4) + 0 + t63 * t55; -t55 * t61 - t56 * t64, -t55 * t62 + t56 * t65, -t56 * t59, -pkin(2) * cos(t58) - cos(qJ(1)) * pkin(1) - t55 * qJ(4) + 0 - t63 * t56; 0, 0, 0, 1;];
	Tc_mdh = t1;
end