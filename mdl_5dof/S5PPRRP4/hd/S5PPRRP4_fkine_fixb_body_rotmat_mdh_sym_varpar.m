% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRRP4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:54
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:44
	% EndTime: 2020-11-04 19:54:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:44
	% EndTime: 2020-11-04 19:54:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t37 = cos(pkin(7));
	t36 = sin(pkin(7));
	t1 = [t37, -t36, 0, 0; t36, t37, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:44
	% EndTime: 2020-11-04 19:54:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t39 = cos(pkin(7));
	t38 = sin(pkin(7));
	t1 = [t39, 0, t38, t39 * pkin(1) + t38 * qJ(2) + 0; t38, 0, -t39, t38 * pkin(1) - t39 * qJ(2) + 0; 0, 1, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:44
	% EndTime: 2020-11-04 19:54:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t46 = pkin(1) + pkin(2);
	t45 = cos(qJ(3));
	t44 = sin(qJ(3));
	t43 = cos(pkin(7));
	t42 = sin(pkin(7));
	t41 = t42 * t45 - t43 * t44;
	t40 = -t42 * t44 - t43 * t45;
	t1 = [-t40, t41, 0, t42 * qJ(2) + t46 * t43 + 0; t41, t40, 0, -t43 * qJ(2) + t46 * t42 + 0; 0, 0, -1, -pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:44
	% EndTime: 2020-11-04 19:54:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (25->19), mult. (32->20), div. (0->0), fcn. (46->6), ass. (0->12)
	t57 = pkin(1) + pkin(2);
	t56 = cos(qJ(3));
	t55 = cos(qJ(4));
	t54 = sin(qJ(3));
	t53 = sin(qJ(4));
	t52 = cos(pkin(7));
	t51 = sin(pkin(7));
	t50 = t52 * pkin(3) - pkin(6) * t51;
	t49 = pkin(3) * t51 + t52 * pkin(6);
	t48 = t51 * t54 + t52 * t56;
	t47 = -t51 * t56 + t52 * t54;
	t1 = [t48 * t55, -t48 * t53, t47, t51 * qJ(2) + t49 * t54 + t50 * t56 + t57 * t52 + 0; -t47 * t55, t47 * t53, t48, -t52 * qJ(2) + t49 * t56 - t50 * t54 + t57 * t51 + 0; -t53, -t55, 0, -pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:54:44
	% EndTime: 2020-11-04 19:54:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->20), mult. (35->18), div. (0->0), fcn. (53->6), ass. (0->12)
	t69 = cos(qJ(3));
	t68 = pkin(1) + pkin(2);
	t67 = cos(qJ(4));
	t66 = sin(qJ(3));
	t65 = sin(qJ(4));
	t64 = qJ(5) + pkin(6);
	t63 = cos(pkin(7));
	t62 = sin(pkin(7));
	t61 = t67 * pkin(4) + pkin(3);
	t59 = t62 * t66 + t63 * t69;
	t58 = -t62 * t69 + t63 * t66;
	t1 = [t59 * t67, -t59 * t65, t58, t62 * qJ(2) + t58 * t64 + t59 * t61 + t68 * t63 + 0; -t58 * t67, t58 * t65, t59, -t63 * qJ(2) - t58 * t61 + t59 * t64 + t68 * t62 + 0; -t65, -t67, 0, -t65 * pkin(4) - pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end