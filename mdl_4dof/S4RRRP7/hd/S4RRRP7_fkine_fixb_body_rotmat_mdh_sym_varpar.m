% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRRP7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:50
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4RRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:50:06
	% EndTime: 2020-11-04 19:50:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:50:06
	% EndTime: 2020-11-04 19:50:06
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t35 = cos(qJ(1));
	t34 = sin(qJ(1));
	t1 = [t35, -t34, 0, 0; t34, t35, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:50:06
	% EndTime: 2020-11-04 19:50:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t39 = cos(qJ(1));
	t38 = cos(qJ(2));
	t37 = sin(qJ(1));
	t36 = sin(qJ(2));
	t1 = [t39 * t38, -t39 * t36, t37, t39 * pkin(1) + t37 * pkin(5) + 0; t37 * t38, -t37 * t36, -t39, t37 * pkin(1) - t39 * pkin(5) + 0; t36, t38, 0, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:50:06
	% EndTime: 2020-11-04 19:50:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t43 = sin(qJ(1));
	t45 = cos(qJ(2));
	t49 = t43 * t45;
	t41 = sin(qJ(3));
	t46 = cos(qJ(1));
	t48 = t46 * t41;
	t44 = cos(qJ(3));
	t47 = t46 * t44;
	t42 = sin(qJ(2));
	t40 = t45 * pkin(2) + t42 * pkin(6) + pkin(1);
	t1 = [t43 * t41 + t45 * t47, t43 * t44 - t45 * t48, t46 * t42, t43 * pkin(5) + t40 * t46 + 0; t44 * t49 - t48, -t41 * t49 - t47, t43 * t42, -t46 * pkin(5) + t40 * t43 + 0; t42 * t44, -t42 * t41, -t45, t42 * pkin(2) - t45 * pkin(6) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:50:06
	% EndTime: 2020-11-04 19:50:06
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (26->18), mult. (36->25), div. (0->0), fcn. (49->6), ass. (0->13)
	t55 = sin(qJ(1));
	t57 = cos(qJ(2));
	t61 = t55 * t57;
	t53 = sin(qJ(3));
	t58 = cos(qJ(1));
	t60 = t58 * t53;
	t56 = cos(qJ(3));
	t59 = t58 * t56;
	t54 = sin(qJ(2));
	t52 = -t53 * pkin(3) + qJ(4) * t56 - pkin(5);
	t51 = t56 * pkin(3) + t53 * qJ(4) + pkin(2);
	t50 = t54 * pkin(6) + t51 * t57 + pkin(1);
	t1 = [t55 * t53 + t57 * t59, t58 * t54, -t55 * t56 + t57 * t60, t50 * t58 - t52 * t55 + 0; t56 * t61 - t60, t55 * t54, t53 * t61 + t59, t50 * t55 + t52 * t58 + 0; t54 * t56, -t57, t54 * t53, -t57 * pkin(6) + t51 * t54 + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end