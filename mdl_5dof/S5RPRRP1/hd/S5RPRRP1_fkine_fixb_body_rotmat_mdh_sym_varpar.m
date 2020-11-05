% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP1 (for one body)
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
% Datum: 2020-11-04 20:22
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:37
	% EndTime: 2020-11-04 20:22:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:37
	% EndTime: 2020-11-04 20:22:37
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t34 = cos(qJ(1));
	t33 = sin(qJ(1));
	t1 = [t34, -t33, 0, 0; t33, t34, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:37
	% EndTime: 2020-11-04 20:22:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t36 = cos(qJ(1));
	t35 = sin(qJ(1));
	t1 = [0, -t36, t35, t36 * pkin(1) + t35 * qJ(2) + 0; 0, -t35, -t36, t35 * pkin(1) - t36 * qJ(2) + 0; 1, 0, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:37
	% EndTime: 2020-11-04 20:22:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t41 = pkin(1) + pkin(6);
	t40 = cos(qJ(1));
	t39 = cos(qJ(3));
	t38 = sin(qJ(1));
	t37 = sin(qJ(3));
	t1 = [t38 * t37, t38 * t39, t40, t38 * qJ(2) + t41 * t40 + 0; -t40 * t37, -t40 * t39, t38, -t40 * qJ(2) + t41 * t38 + 0; t39, -t37, 0, pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:37
	% EndTime: 2020-11-04 20:22:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->13), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t48 = cos(qJ(1));
	t47 = sin(qJ(1));
	t46 = qJ(3) + qJ(4);
	t45 = pkin(1) + pkin(6) + pkin(7);
	t44 = cos(t46);
	t43 = sin(t46);
	t42 = sin(qJ(3)) * pkin(3) + qJ(2);
	t1 = [t47 * t43, t47 * t44, t48, t42 * t47 + t45 * t48 + 0; -t48 * t43, -t48 * t44, t47, -t42 * t48 + t45 * t47 + 0; t44, -t43, 0, cos(qJ(3)) * pkin(3) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:22:37
	% EndTime: 2020-11-04 20:22:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->15), mult. (16->12), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = qJ(3) + qJ(4);
	t50 = sin(t53);
	t56 = pkin(4) * t50 + sin(qJ(3)) * pkin(3) + qJ(2);
	t55 = cos(qJ(1));
	t54 = sin(qJ(1));
	t52 = qJ(5) + pkin(1) + pkin(6) + pkin(7);
	t51 = cos(t53);
	t1 = [t54 * t50, t54 * t51, t55, t52 * t55 + t56 * t54 + 0; -t55 * t50, -t55 * t51, t54, t52 * t54 - t56 * t55 + 0; t51, -t50, 0, pkin(4) * t51 + cos(qJ(3)) * pkin(3) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end