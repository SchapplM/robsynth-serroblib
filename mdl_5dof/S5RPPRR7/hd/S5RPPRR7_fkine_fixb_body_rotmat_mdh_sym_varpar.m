% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:15
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPRR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:19
	% EndTime: 2020-11-04 20:15:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:19
	% EndTime: 2020-11-04 20:15:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t38 = cos(qJ(1));
	t37 = sin(qJ(1));
	t1 = [t38, -t37, 0, 0; t37, t38, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:19
	% EndTime: 2020-11-04 20:15:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t41 = qJ(1) + pkin(8);
	t40 = cos(t41);
	t39 = sin(t41);
	t1 = [t40, -t39, 0, cos(qJ(1)) * pkin(1) + 0; t39, t40, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:19
	% EndTime: 2020-11-04 20:15:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (6->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t44 = qJ(1) + pkin(8);
	t43 = cos(t44);
	t42 = sin(t44);
	t1 = [0, -t43, t42, t43 * pkin(2) + t42 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; 0, -t42, -t43, t42 * pkin(2) - t43 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; 1, 0, 0, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:19
	% EndTime: 2020-11-04 20:15:19
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (24->14), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->7)
	t50 = pkin(2) + pkin(6);
	t49 = cos(qJ(4));
	t48 = sin(qJ(4));
	t47 = qJ(1) + pkin(8);
	t46 = cos(t47);
	t45 = sin(t47);
	t1 = [t45 * t48, t45 * t49, t46, t50 * t46 + t45 * qJ(3) + cos(qJ(1)) * pkin(1) + 0; -t46 * t48, -t46 * t49, t45, t50 * t45 - t46 * qJ(3) + sin(qJ(1)) * pkin(1) + 0; t49, -t48, 0, pkin(3) + qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:15:19
	% EndTime: 2020-11-04 20:15:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (41->21), mult. (32->24), div. (0->0), fcn. (45->8), ass. (0->12)
	t54 = sin(qJ(5));
	t55 = sin(qJ(4));
	t61 = t54 * t55;
	t56 = cos(qJ(5));
	t60 = t55 * t56;
	t57 = cos(qJ(4));
	t59 = pkin(4) * t55 - pkin(7) * t57 + qJ(3);
	t58 = pkin(2) + pkin(6);
	t53 = qJ(1) + pkin(8);
	t52 = cos(t53);
	t51 = sin(t53);
	t1 = [t51 * t60 + t52 * t54, -t51 * t61 + t52 * t56, -t51 * t57, cos(qJ(1)) * pkin(1) + t58 * t52 + 0 + t59 * t51; t51 * t54 - t52 * t60, t51 * t56 + t52 * t61, t52 * t57, sin(qJ(1)) * pkin(1) + t58 * t51 + 0 - t59 * t52; t57 * t56, -t57 * t54, t55, t57 * pkin(4) + t55 * pkin(7) + pkin(3) + pkin(5) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end