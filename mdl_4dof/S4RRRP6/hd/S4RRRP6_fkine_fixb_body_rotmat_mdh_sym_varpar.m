% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRRP6 (for one body)
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
% Datum: 2020-11-04 19:49
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4RRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:49:52
	% EndTime: 2020-11-04 19:49:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:49:52
	% EndTime: 2020-11-04 19:49:52
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t33 = cos(qJ(1));
	t32 = sin(qJ(1));
	t1 = [t33, -t32, 0, 0; t32, t33, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:49:52
	% EndTime: 2020-11-04 19:49:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t37 = cos(qJ(1));
	t36 = cos(qJ(2));
	t35 = sin(qJ(1));
	t34 = sin(qJ(2));
	t1 = [t37 * t36, -t37 * t34, t35, t37 * pkin(1) + t35 * pkin(5) + 0; t35 * t36, -t35 * t34, -t37, t35 * pkin(1) - t37 * pkin(5) + 0; t34, t36, 0, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:49:52
	% EndTime: 2020-11-04 19:49:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t41 = sin(qJ(1));
	t43 = cos(qJ(2));
	t47 = t41 * t43;
	t39 = sin(qJ(3));
	t44 = cos(qJ(1));
	t46 = t44 * t39;
	t42 = cos(qJ(3));
	t45 = t44 * t42;
	t40 = sin(qJ(2));
	t38 = t43 * pkin(2) + t40 * pkin(6) + pkin(1);
	t1 = [t41 * t39 + t43 * t45, t41 * t42 - t43 * t46, t44 * t40, t41 * pkin(5) + t38 * t44 + 0; t42 * t47 - t46, -t39 * t47 - t45, t41 * t40, -t44 * pkin(5) + t38 * t41 + 0; t40 * t42, -t40 * t39, -t43, t40 * pkin(2) - t43 * pkin(6) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:49:52
	% EndTime: 2020-11-04 19:49:52
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (25->18), mult. (31->23), div. (0->0), fcn. (44->6), ass. (0->14)
	t54 = sin(qJ(1));
	t56 = cos(qJ(2));
	t60 = t54 * t56;
	t52 = sin(qJ(3));
	t57 = cos(qJ(1));
	t59 = t57 * t52;
	t55 = cos(qJ(3));
	t58 = t57 * t55;
	t53 = sin(qJ(2));
	t51 = qJ(4) + pkin(6);
	t50 = t55 * pkin(3) + pkin(2);
	t49 = t52 * pkin(3) + pkin(5);
	t48 = t50 * t56 + t51 * t53 + pkin(1);
	t1 = [t54 * t52 + t56 * t58, t54 * t55 - t56 * t59, t57 * t53, t48 * t57 + t49 * t54 + 0; t55 * t60 - t59, -t52 * t60 - t58, t54 * t53, t48 * t54 - t49 * t57 + 0; t53 * t55, -t53 * t52, -t56, t53 * t50 - t56 * t51 + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end