% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRR8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:27
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRRR8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:27:31
	% EndTime: 2020-11-04 20:27:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:27:31
	% EndTime: 2020-11-04 20:27:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t36 = cos(qJ(1));
	t35 = sin(qJ(1));
	t1 = [t36, -t35, 0, 0; t35, t36, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:27:31
	% EndTime: 2020-11-04 20:27:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t38 = cos(qJ(1));
	t37 = sin(qJ(1));
	t1 = [t38, 0, t37, t38 * pkin(1) + t37 * qJ(2) + 0; t37, 0, -t38, t37 * pkin(1) - t38 * qJ(2) + 0; 0, 1, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:27:31
	% EndTime: 2020-11-04 20:27:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t45 = pkin(1) + pkin(2);
	t44 = cos(qJ(1));
	t43 = cos(qJ(3));
	t42 = sin(qJ(1));
	t41 = sin(qJ(3));
	t40 = -t44 * t41 + t42 * t43;
	t39 = -t42 * t41 - t44 * t43;
	t1 = [-t39, t40, 0, t42 * qJ(2) + t45 * t44 + 0; t40, t39, 0, -t44 * qJ(2) + t45 * t42 + 0; 0, 0, -1, -pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:27:31
	% EndTime: 2020-11-04 20:27:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (25->17), mult. (28->16), div. (0->0), fcn. (42->6), ass. (0->11)
	t55 = cos(qJ(1));
	t54 = cos(qJ(3));
	t53 = cos(qJ(4));
	t52 = sin(qJ(1));
	t51 = sin(qJ(3));
	t50 = sin(qJ(4));
	t49 = -pkin(3) * t51 + pkin(7) * t54 - qJ(2);
	t48 = pkin(3) * t54 + pkin(7) * t51 + pkin(1) + pkin(2);
	t47 = t51 * t52 + t54 * t55;
	t46 = t51 * t55 - t52 * t54;
	t1 = [t47 * t53, -t47 * t50, t46, t48 * t55 - t49 * t52 + 0; -t46 * t53, t46 * t50, t47, t48 * t52 + t49 * t55 + 0; -t50, -t53, 0, -pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:27:31
	% EndTime: 2020-11-04 20:27:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (40->21), mult. (33->18), div. (0->0), fcn. (47->8), ass. (0->14)
	t59 = cos(qJ(4)) * pkin(4) + pkin(3);
	t63 = sin(qJ(3));
	t65 = cos(qJ(3));
	t67 = pkin(8) + pkin(7);
	t68 = t59 * t63 - t67 * t65 + qJ(2);
	t66 = cos(qJ(1));
	t64 = sin(qJ(1));
	t62 = qJ(4) + qJ(5);
	t61 = cos(t62);
	t60 = sin(t62);
	t58 = t64 * t63 + t66 * t65;
	t57 = t66 * t63 - t64 * t65;
	t56 = t59 * t65 + t67 * t63 + pkin(1) + pkin(2);
	t1 = [t58 * t61, -t58 * t60, t57, t56 * t66 + t68 * t64 + 0; -t57 * t61, t57 * t60, t58, t56 * t64 - t68 * t66 + 0; -t60, -t61, 0, -sin(qJ(4)) * pkin(4) - pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end