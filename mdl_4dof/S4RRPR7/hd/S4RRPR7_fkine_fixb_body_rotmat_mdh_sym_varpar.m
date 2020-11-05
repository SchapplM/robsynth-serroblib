% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRPR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:47
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4RRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:47:45
	% EndTime: 2020-11-04 19:47:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:47:45
	% EndTime: 2020-11-04 19:47:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t37 = cos(qJ(1));
	t36 = sin(qJ(1));
	t1 = [t37, -t36, 0, 0; t36, t37, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:47:45
	% EndTime: 2020-11-04 19:47:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t41 = cos(qJ(1));
	t40 = cos(qJ(2));
	t39 = sin(qJ(1));
	t38 = sin(qJ(2));
	t1 = [t41 * t40, -t41 * t38, t39, t41 * pkin(1) + t39 * pkin(5) + 0; t39 * t40, -t39 * t38, -t41, t39 * pkin(1) - t41 * pkin(5) + 0; t38, t40, 0, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:47:45
	% EndTime: 2020-11-04 19:47:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t48 = cos(qJ(1));
	t47 = sin(qJ(1));
	t46 = -qJ(3) - pkin(5);
	t45 = qJ(2) + pkin(7);
	t44 = cos(t45);
	t43 = sin(t45);
	t42 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t48 * t44, -t48 * t43, t47, t48 * t42 - t46 * t47 + 0; t47 * t44, -t47 * t43, -t48, t47 * t42 + t48 * t46 + 0; t43, t44, 0, sin(qJ(2)) * pkin(2) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:47:45
	% EndTime: 2020-11-04 19:47:45
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (37->21), mult. (35->25), div. (0->0), fcn. (48->10), ass. (0->17)
	t56 = sin(qJ(4));
	t58 = sin(qJ(1));
	t64 = t58 * t56;
	t59 = cos(qJ(4));
	t63 = t58 * t59;
	t60 = cos(qJ(1));
	t62 = t60 * t56;
	t61 = t60 * t59;
	t57 = sin(qJ(2));
	t55 = -qJ(3) - pkin(5);
	t54 = cos(pkin(7));
	t53 = sin(pkin(7));
	t52 = qJ(2) + pkin(7);
	t51 = cos(t52);
	t50 = sin(t52);
	t49 = (pkin(3) * t54 + pkin(6) * t53 + pkin(2)) * cos(qJ(2)) + (-pkin(3) * t53 + pkin(6) * t54) * t57 + pkin(1);
	t1 = [t51 * t61 + t64, -t51 * t62 + t63, t60 * t50, t49 * t60 - t55 * t58 + 0; t51 * t63 - t62, -t51 * t64 - t61, t58 * t50, t49 * t58 + t55 * t60 + 0; t50 * t59, -t50 * t56, -t51, pkin(2) * t57 + pkin(3) * t50 - pkin(6) * t51 + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end