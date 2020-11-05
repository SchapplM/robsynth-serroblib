% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:51
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4RRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:51:24
	% EndTime: 2020-11-04 19:51:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:51:24
	% EndTime: 2020-11-04 19:51:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t32 = cos(qJ(1));
	t31 = sin(qJ(1));
	t1 = [t32, -t31, 0, 0; t31, t32, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:51:24
	% EndTime: 2020-11-04 19:51:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t36 = cos(qJ(1));
	t35 = cos(qJ(2));
	t34 = sin(qJ(1));
	t33 = sin(qJ(2));
	t1 = [t36 * t35, -t36 * t33, t34, t36 * pkin(1) + t34 * pkin(5) + 0; t34 * t35, -t34 * t33, -t36, t34 * pkin(1) - t36 * pkin(5) + 0; t33, t35, 0, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:51:24
	% EndTime: 2020-11-04 19:51:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t40 = sin(qJ(1));
	t42 = cos(qJ(2));
	t46 = t40 * t42;
	t38 = sin(qJ(3));
	t43 = cos(qJ(1));
	t45 = t43 * t38;
	t41 = cos(qJ(3));
	t44 = t43 * t41;
	t39 = sin(qJ(2));
	t37 = t42 * pkin(2) + t39 * pkin(6) + pkin(1);
	t1 = [t40 * t38 + t42 * t44, t40 * t41 - t42 * t45, t43 * t39, t40 * pkin(5) + t37 * t43 + 0; t41 * t46 - t45, -t38 * t46 - t44, t40 * t39, -t43 * pkin(5) + t37 * t40 + 0; t39 * t41, -t39 * t38, -t42, t39 * pkin(2) - t42 * pkin(6) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:51:24
	% EndTime: 2020-11-04 19:51:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (31->23), div. (0->0), fcn. (44->8), ass. (0->15)
	t54 = sin(qJ(1));
	t55 = cos(qJ(2));
	t60 = t54 * t55;
	t52 = qJ(3) + qJ(4);
	t50 = sin(t52);
	t56 = cos(qJ(1));
	t59 = t56 * t50;
	t51 = cos(t52);
	t58 = t56 * t51;
	t57 = pkin(7) + pkin(6);
	t53 = sin(qJ(2));
	t49 = cos(qJ(3)) * pkin(3) + pkin(2);
	t48 = sin(qJ(3)) * pkin(3) + pkin(5);
	t47 = t49 * t55 + t57 * t53 + pkin(1);
	t1 = [t54 * t50 + t55 * t58, t54 * t51 - t55 * t59, t56 * t53, t47 * t56 + t48 * t54 + 0; t51 * t60 - t59, -t50 * t60 - t58, t54 * t53, t47 * t54 - t48 * t56 + 0; t53 * t51, -t53 * t50, -t55, t53 * t49 - t55 * t57 + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end