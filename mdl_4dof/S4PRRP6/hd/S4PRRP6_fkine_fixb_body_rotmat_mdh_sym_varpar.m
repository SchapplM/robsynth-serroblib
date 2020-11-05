% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRRP6 (for one body)
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

function Tc_mdh = S4PRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:37:25
	% EndTime: 2020-11-04 19:37:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:37:25
	% EndTime: 2020-11-04 19:37:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t39 = cos(pkin(6));
	t38 = sin(pkin(6));
	t1 = [t39, -t38, 0, 0; t38, t39, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:37:25
	% EndTime: 2020-11-04 19:37:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t43 = cos(qJ(2));
	t42 = sin(qJ(2));
	t41 = cos(pkin(6));
	t40 = sin(pkin(6));
	t1 = [t41 * t43, -t41 * t42, t40, t41 * pkin(1) + t40 * pkin(4) + 0; t40 * t43, -t40 * t42, -t41, t40 * pkin(1) - t41 * pkin(4) + 0; t42, t43, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:37:25
	% EndTime: 2020-11-04 19:37:25
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (30->22), div. (0->0), fcn. (43->6), ass. (0->10)
	t46 = sin(qJ(3));
	t49 = cos(qJ(2));
	t52 = t46 * t49;
	t48 = cos(qJ(3));
	t51 = t48 * t49;
	t47 = sin(qJ(2));
	t50 = pkin(2) * t49 + pkin(5) * t47 + pkin(1);
	t45 = cos(pkin(6));
	t44 = sin(pkin(6));
	t1 = [t44 * t46 + t45 * t51, t44 * t48 - t45 * t52, t45 * t47, t44 * pkin(4) + t50 * t45 + 0; t44 * t51 - t45 * t46, -t44 * t52 - t45 * t48, t44 * t47, -t45 * pkin(4) + t50 * t44 + 0; t47 * t48, -t47 * t46, -t49, t47 * pkin(2) - t49 * pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:37:25
	% EndTime: 2020-11-04 19:37:26
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (26->18), mult. (44->26), div. (0->0), fcn. (57->6), ass. (0->12)
	t56 = sin(qJ(3));
	t59 = cos(qJ(2));
	t63 = t56 * t59;
	t58 = cos(qJ(3));
	t62 = t58 * t59;
	t53 = pkin(3) * t58 + qJ(4) * t56 + pkin(2);
	t57 = sin(qJ(2));
	t61 = pkin(5) * t57 + t53 * t59 + pkin(1);
	t60 = pkin(3) * t56 - qJ(4) * t58 + pkin(4);
	t55 = cos(pkin(6));
	t54 = sin(pkin(6));
	t1 = [t54 * t56 + t55 * t62, t55 * t57, -t54 * t58 + t55 * t63, t54 * t60 + t55 * t61 + 0; t54 * t62 - t55 * t56, t54 * t57, t54 * t63 + t55 * t58, t54 * t61 - t55 * t60 + 0; t57 * t58, -t59, t57 * t56, -pkin(5) * t59 + t53 * t57 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end