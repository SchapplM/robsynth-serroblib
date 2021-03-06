% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRPR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:35
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4PRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:35:52
	% EndTime: 2020-11-04 19:35:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:35:52
	% EndTime: 2020-11-04 19:35:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t35 = cos(pkin(6));
	t34 = sin(pkin(6));
	t1 = [t35, -t34, 0, 0; t34, t35, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:35:52
	% EndTime: 2020-11-04 19:35:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t39 = cos(qJ(2));
	t38 = sin(qJ(2));
	t37 = cos(pkin(6));
	t36 = sin(pkin(6));
	t1 = [t37 * t39, -t37 * t38, t36, t37 * pkin(1) + t36 * pkin(4) + 0; t36 * t39, -t36 * t38, -t37, t36 * pkin(1) - t37 * pkin(4) + 0; t38, t39, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:35:52
	% EndTime: 2020-11-04 19:35:52
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (16->14), mult. (18->12), div. (0->0), fcn. (26->4), ass. (0->6)
	t42 = sin(qJ(2));
	t43 = cos(qJ(2));
	t44 = pkin(2) * t43 + qJ(3) * t42 + pkin(1);
	t41 = cos(pkin(6));
	t40 = sin(pkin(6));
	t1 = [t40, -t41 * t43, t41 * t42, t40 * pkin(4) + t44 * t41 + 0; -t41, -t40 * t43, t40 * t42, -t41 * pkin(4) + t44 * t40 + 0; 0, -t42, -t43, t42 * pkin(2) - t43 * qJ(3) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:35:52
	% EndTime: 2020-11-04 19:35:52
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->17), mult. (30->22), div. (0->0), fcn. (43->6), ass. (0->12)
	t47 = sin(qJ(4));
	t48 = sin(qJ(2));
	t55 = t47 * t48;
	t49 = cos(qJ(4));
	t54 = t48 * t49;
	t50 = cos(qJ(2));
	t52 = pkin(2) + pkin(5);
	t53 = qJ(3) * t48 + t50 * t52 + pkin(1);
	t51 = pkin(3) + pkin(4);
	t46 = cos(pkin(6));
	t45 = sin(pkin(6));
	t1 = [t45 * t49 + t46 * t55, -t45 * t47 + t46 * t54, t46 * t50, t51 * t45 + t53 * t46 + 0; t45 * t55 - t46 * t49, t45 * t54 + t46 * t47, t45 * t50, t53 * t45 - t51 * t46 + 0; -t50 * t47, -t50 * t49, t48, -t50 * qJ(3) + t52 * t48 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end