% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:24
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:24:31
	% EndTime: 2020-11-04 20:24:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:24:31
	% EndTime: 2020-11-04 20:24:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t37 = cos(qJ(1));
	t36 = sin(qJ(1));
	t1 = [t37, -t36, 0, 0; t36, t37, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:24:31
	% EndTime: 2020-11-04 20:24:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t39 = cos(qJ(1));
	t38 = sin(qJ(1));
	t1 = [t39, 0, t38, t39 * pkin(1) + t38 * qJ(2) + 0; t38, 0, -t39, t38 * pkin(1) - t39 * qJ(2) + 0; 0, 1, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:24:31
	% EndTime: 2020-11-04 20:24:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t46 = pkin(1) + pkin(2);
	t45 = cos(qJ(1));
	t44 = cos(qJ(3));
	t43 = sin(qJ(1));
	t42 = sin(qJ(3));
	t41 = -t45 * t42 + t43 * t44;
	t40 = -t43 * t42 - t45 * t44;
	t1 = [-t40, t41, 0, t43 * qJ(2) + t46 * t45 + 0; t41, t40, 0, -t45 * qJ(2) + t46 * t43 + 0; 0, 0, -1, -pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:24:31
	% EndTime: 2020-11-04 20:24:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (25->17), mult. (28->16), div. (0->0), fcn. (42->6), ass. (0->11)
	t56 = cos(qJ(1));
	t55 = cos(qJ(3));
	t54 = cos(qJ(4));
	t53 = sin(qJ(1));
	t52 = sin(qJ(3));
	t51 = sin(qJ(4));
	t50 = -t52 * pkin(3) + t55 * pkin(7) - qJ(2);
	t49 = t55 * pkin(3) + t52 * pkin(7) + pkin(1) + pkin(2);
	t48 = t53 * t52 + t56 * t55;
	t47 = t56 * t52 - t53 * t55;
	t1 = [t48 * t54, -t48 * t51, t47, t49 * t56 - t50 * t53 + 0; -t47 * t54, t47 * t51, t48, t49 * t53 + t50 * t56 + 0; -t51, -t54, 0, -pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:24:31
	% EndTime: 2020-11-04 20:24:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (33->20), mult. (38->20), div. (0->0), fcn. (52->6), ass. (0->12)
	t62 = sin(qJ(4));
	t65 = cos(qJ(4));
	t60 = pkin(4) * t65 + qJ(5) * t62 + pkin(3);
	t63 = sin(qJ(3));
	t66 = cos(qJ(3));
	t68 = t66 * pkin(7) - t60 * t63 - qJ(2);
	t67 = cos(qJ(1));
	t64 = sin(qJ(1));
	t59 = t64 * t63 + t67 * t66;
	t58 = t67 * t63 - t64 * t66;
	t57 = t63 * pkin(7) + t60 * t66 + pkin(1) + pkin(2);
	t1 = [t59 * t65, t58, t59 * t62, t57 * t67 - t68 * t64 + 0; -t58 * t65, t59, -t58 * t62, t57 * t64 + t68 * t67 + 0; -t62, 0, t65, -t62 * pkin(4) + t65 * qJ(5) + pkin(5) - pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end