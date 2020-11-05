% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPRR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:43
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4RPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:43:47
	% EndTime: 2020-11-04 19:43:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:43:47
	% EndTime: 2020-11-04 19:43:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t35 = cos(qJ(1));
	t34 = sin(qJ(1));
	t1 = [t35, -t34, 0, 0; t34, t35, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:43:47
	% EndTime: 2020-11-04 19:43:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t38 = qJ(1) + pkin(7);
	t37 = cos(t38);
	t36 = sin(t38);
	t1 = [t37, -t36, 0, cos(qJ(1)) * pkin(1) + 0; t36, t37, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:43:47
	% EndTime: 2020-11-04 19:43:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t43 = cos(qJ(3));
	t42 = sin(qJ(3));
	t41 = qJ(1) + pkin(7);
	t40 = cos(t41);
	t39 = sin(t41);
	t1 = [t40 * t43, -t40 * t42, t39, t40 * pkin(2) + t39 * pkin(5) + cos(qJ(1)) * pkin(1) + 0; t39 * t43, -t39 * t42, -t40, t39 * pkin(2) - t40 * pkin(5) + sin(qJ(1)) * pkin(1) + 0; t42, t43, 0, qJ(2) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:43:47
	% EndTime: 2020-11-04 19:43:47
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (38->19), mult. (32->24), div. (0->0), fcn. (45->8), ass. (0->11)
	t47 = sin(qJ(4));
	t50 = cos(qJ(3));
	t53 = t47 * t50;
	t49 = cos(qJ(4));
	t52 = t49 * t50;
	t48 = sin(qJ(3));
	t51 = pkin(3) * t50 + pkin(6) * t48 + pkin(2);
	t46 = qJ(1) + pkin(7);
	t45 = cos(t46);
	t44 = sin(t46);
	t1 = [t44 * t47 + t45 * t52, t44 * t49 - t45 * t53, t45 * t48, cos(qJ(1)) * pkin(1) + t44 * pkin(5) + 0 + t51 * t45; t44 * t52 - t45 * t47, -t44 * t53 - t45 * t49, t44 * t48, sin(qJ(1)) * pkin(1) - t45 * pkin(5) + 0 + t51 * t44; t48 * t49, -t48 * t47, -t50, pkin(3) * t48 - pkin(6) * t50 + pkin(4) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end