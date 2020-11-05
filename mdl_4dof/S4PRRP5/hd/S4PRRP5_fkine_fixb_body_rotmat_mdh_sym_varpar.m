% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRRP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:37
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4PRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:37:11
	% EndTime: 2020-11-04 19:37:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:37:11
	% EndTime: 2020-11-04 19:37:11
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t36 = cos(pkin(6));
	t35 = sin(pkin(6));
	t1 = [t36, -t35, 0, 0; t35, t36, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:37:11
	% EndTime: 2020-11-04 19:37:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t40 = cos(qJ(2));
	t39 = sin(qJ(2));
	t38 = cos(pkin(6));
	t37 = sin(pkin(6));
	t1 = [t38 * t40, -t38 * t39, t37, t38 * pkin(1) + t37 * pkin(4) + 0; t37 * t40, -t37 * t39, -t38, t37 * pkin(1) - t38 * pkin(4) + 0; t39, t40, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:37:11
	% EndTime: 2020-11-04 19:37:11
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (30->22), div. (0->0), fcn. (43->6), ass. (0->10)
	t43 = sin(qJ(3));
	t46 = cos(qJ(2));
	t49 = t43 * t46;
	t45 = cos(qJ(3));
	t48 = t45 * t46;
	t44 = sin(qJ(2));
	t47 = pkin(2) * t46 + pkin(5) * t44 + pkin(1);
	t42 = cos(pkin(6));
	t41 = sin(pkin(6));
	t1 = [t41 * t43 + t42 * t48, t41 * t45 - t42 * t49, t42 * t44, t41 * pkin(4) + t47 * t42 + 0; t41 * t48 - t42 * t43, -t41 * t49 - t42 * t45, t41 * t44, -t42 * pkin(4) + t47 * t41 + 0; t44 * t45, -t44 * t43, -t46, t44 * pkin(2) - t46 * pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:37:11
	% EndTime: 2020-11-04 19:37:11
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (25->18), mult. (37->24), div. (0->0), fcn. (50->6), ass. (0->13)
	t54 = sin(qJ(3));
	t57 = cos(qJ(2));
	t61 = t54 * t57;
	t56 = cos(qJ(3));
	t60 = t56 * t57;
	t59 = pkin(3) * t54 + pkin(4);
	t50 = t56 * pkin(3) + pkin(2);
	t53 = qJ(4) + pkin(5);
	t55 = sin(qJ(2));
	t58 = t50 * t57 + t53 * t55 + pkin(1);
	t52 = cos(pkin(6));
	t51 = sin(pkin(6));
	t1 = [t51 * t54 + t52 * t60, t51 * t56 - t52 * t61, t52 * t55, t59 * t51 + t58 * t52 + 0; t51 * t60 - t52 * t54, -t51 * t61 - t52 * t56, t51 * t55, t58 * t51 - t59 * t52 + 0; t55 * t56, -t55 * t54, -t57, t55 * t50 - t57 * t53 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end