% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRRR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:39
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4PRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:39:01
	% EndTime: 2020-11-04 19:39:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:39:01
	% EndTime: 2020-11-04 19:39:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t37 = cos(pkin(7));
	t36 = sin(pkin(7));
	t1 = [t37, -t36, 0, 0; t36, t37, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:39:01
	% EndTime: 2020-11-04 19:39:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t41 = cos(qJ(2));
	t40 = sin(qJ(2));
	t39 = cos(pkin(7));
	t38 = sin(pkin(7));
	t1 = [t39 * t41, -t39 * t40, t38, t39 * pkin(1) + t38 * pkin(4) + 0; t38 * t41, -t38 * t40, -t39, t38 * pkin(1) - t39 * pkin(4) + 0; t40, t41, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:39:01
	% EndTime: 2020-11-04 19:39:01
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (30->22), div. (0->0), fcn. (43->6), ass. (0->10)
	t44 = sin(qJ(3));
	t47 = cos(qJ(2));
	t50 = t44 * t47;
	t46 = cos(qJ(3));
	t49 = t46 * t47;
	t45 = sin(qJ(2));
	t48 = pkin(2) * t47 + pkin(5) * t45 + pkin(1);
	t43 = cos(pkin(7));
	t42 = sin(pkin(7));
	t1 = [t42 * t44 + t43 * t49, t42 * t46 - t43 * t50, t43 * t45, t42 * pkin(4) + t48 * t43 + 0; t42 * t49 - t43 * t44, -t42 * t50 - t43 * t46, t42 * t45, -t43 * pkin(4) + t48 * t42 + 0; t45 * t46, -t45 * t44, -t47, t45 * pkin(2) - t47 * pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:39:01
	% EndTime: 2020-11-04 19:39:01
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (35->19), mult. (37->24), div. (0->0), fcn. (50->8), ass. (0->14)
	t55 = sin(pkin(7));
	t59 = cos(qJ(2));
	t64 = t55 * t59;
	t56 = cos(pkin(7));
	t63 = t56 * t59;
	t62 = pkin(3) * sin(qJ(3)) + pkin(4);
	t51 = cos(qJ(3)) * pkin(3) + pkin(2);
	t58 = sin(qJ(2));
	t60 = pkin(6) + pkin(5);
	t61 = t51 * t59 + t58 * t60 + pkin(1);
	t54 = qJ(3) + qJ(4);
	t53 = cos(t54);
	t52 = sin(t54);
	t1 = [t55 * t52 + t53 * t63, -t52 * t63 + t55 * t53, t56 * t58, t62 * t55 + t61 * t56 + 0; -t56 * t52 + t53 * t64, -t52 * t64 - t56 * t53, t55 * t58, t61 * t55 - t62 * t56 + 0; t58 * t53, -t58 * t52, -t59, t58 * t51 - t59 * t60 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end