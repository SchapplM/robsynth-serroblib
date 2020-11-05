% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:38
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4PRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:38:44
	% EndTime: 2020-11-04 19:38:44
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:38:44
	% EndTime: 2020-11-04 19:38:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t36 = cos(pkin(7));
	t35 = sin(pkin(7));
	t1 = [t36, -t35, 0, 0; t35, t36, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:38:44
	% EndTime: 2020-11-04 19:38:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t40 = cos(qJ(2));
	t39 = sin(qJ(2));
	t38 = cos(pkin(7));
	t37 = sin(pkin(7));
	t1 = [t38 * t40, -t38 * t39, t37, t38 * pkin(1) + t37 * pkin(4) + 0; t37 * t40, -t37 * t39, -t38, t37 * pkin(1) - t38 * pkin(4) + 0; t39, t40, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:38:44
	% EndTime: 2020-11-04 19:38:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t47 = -pkin(5) - pkin(4);
	t46 = cos(pkin(7));
	t45 = sin(pkin(7));
	t44 = qJ(2) + qJ(3);
	t43 = cos(t44);
	t42 = sin(t44);
	t41 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t46 * t43, -t46 * t42, t45, t46 * t41 - t45 * t47 + 0; t45 * t43, -t45 * t42, -t46, t45 * t41 + t46 * t47 + 0; t42, t43, 0, sin(qJ(2)) * pkin(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:38:44
	% EndTime: 2020-11-04 19:38:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t52 = sin(pkin(7));
	t54 = sin(qJ(4));
	t61 = t52 * t54;
	t55 = cos(qJ(4));
	t60 = t52 * t55;
	t53 = cos(pkin(7));
	t59 = t53 * t54;
	t58 = t53 * t55;
	t51 = qJ(2) + qJ(3);
	t49 = sin(t51);
	t50 = cos(t51);
	t57 = pkin(3) * t50 + pkin(6) * t49 + cos(qJ(2)) * pkin(2) + pkin(1);
	t56 = -pkin(5) - pkin(4);
	t1 = [t50 * t58 + t61, -t50 * t59 + t60, t53 * t49, -t52 * t56 + t57 * t53 + 0; t50 * t60 - t59, -t50 * t61 - t58, t52 * t49, t57 * t52 + t53 * t56 + 0; t49 * t55, -t49 * t54, -t50, t49 * pkin(3) - t50 * pkin(6) + sin(qJ(2)) * pkin(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end