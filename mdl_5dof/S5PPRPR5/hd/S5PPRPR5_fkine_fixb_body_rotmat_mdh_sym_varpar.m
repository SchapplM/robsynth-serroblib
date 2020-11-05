% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRPR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:53
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:38
	% EndTime: 2020-11-04 19:53:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:38
	% EndTime: 2020-11-04 19:53:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t36 = cos(pkin(7));
	t35 = sin(pkin(7));
	t1 = [t36, -t35, 0, 0; t35, t36, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:38
	% EndTime: 2020-11-04 19:53:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t38 = cos(pkin(7));
	t37 = sin(pkin(7));
	t1 = [t38, 0, t37, t38 * pkin(1) + t37 * qJ(2) + 0; t37, 0, -t38, t37 * pkin(1) - t38 * qJ(2) + 0; 0, 1, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:38
	% EndTime: 2020-11-04 19:53:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t45 = pkin(1) + pkin(2);
	t44 = cos(qJ(3));
	t43 = sin(qJ(3));
	t42 = cos(pkin(7));
	t41 = sin(pkin(7));
	t40 = t41 * t44 - t42 * t43;
	t39 = -t41 * t43 - t42 * t44;
	t1 = [-t39, t40, 0, t41 * qJ(2) + t45 * t42 + 0; t40, t39, 0, -t42 * qJ(2) + t45 * t41 + 0; 0, 0, -1, -pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:38
	% EndTime: 2020-11-04 19:53:38
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (20->16), mult. (24->16), div. (0->0), fcn. (32->4), ass. (0->10)
	t54 = pkin(1) + pkin(2);
	t53 = cos(qJ(3));
	t52 = sin(qJ(3));
	t51 = cos(pkin(7));
	t50 = sin(pkin(7));
	t49 = t51 * pkin(3) - qJ(4) * t50;
	t48 = pkin(3) * t50 + t51 * qJ(4);
	t47 = -t50 * t53 + t51 * t52;
	t46 = -t50 * t52 - t51 * t53;
	t1 = [0, t46, t47, t50 * qJ(2) + t48 * t52 + t49 * t53 + t54 * t51 + 0; 0, t47, -t46, -t51 * qJ(2) + t48 * t53 - t49 * t52 + t54 * t50 + 0; -1, 0, 0, -pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:53:38
	% EndTime: 2020-11-04 19:53:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (28->19), mult. (32->20), div. (0->0), fcn. (46->6), ass. (0->13)
	t66 = pkin(1) + pkin(2);
	t65 = pkin(3) + pkin(6);
	t64 = cos(qJ(3));
	t63 = cos(qJ(5));
	t62 = sin(qJ(3));
	t61 = sin(qJ(5));
	t60 = cos(pkin(7));
	t59 = sin(pkin(7));
	t58 = t60 * qJ(4) + t65 * t59;
	t57 = -qJ(4) * t59 + t60 * t65;
	t56 = t59 * t62 + t60 * t64;
	t55 = -t59 * t64 + t60 * t62;
	t1 = [t55 * t61, t55 * t63, t56, t59 * qJ(2) + t57 * t64 + t58 * t62 + t66 * t60 + 0; t56 * t61, t56 * t63, -t55, -t60 * qJ(2) - t57 * t62 + t58 * t64 + t66 * t59 + 0; -t63, t61, 0, -pkin(4) - pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end