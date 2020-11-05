% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPRR7 (for one body)
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
% Datum: 2020-11-04 19:44
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4RPRR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:44:29
	% EndTime: 2020-11-04 19:44:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:44:29
	% EndTime: 2020-11-04 19:44:29
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t37 = cos(qJ(1));
	t36 = sin(qJ(1));
	t1 = [t37, -t36, 0, 0; t36, t37, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:44:29
	% EndTime: 2020-11-04 19:44:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t41 = cos(qJ(1));
	t40 = sin(qJ(1));
	t39 = cos(pkin(7));
	t38 = sin(pkin(7));
	t1 = [t41 * t39, -t41 * t38, t40, t41 * pkin(1) + t40 * qJ(2) + 0; t40 * t39, -t40 * t38, -t41, t40 * pkin(1) - t41 * qJ(2) + 0; t38, t39, 0, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:44:29
	% EndTime: 2020-11-04 19:44:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t48 = cos(qJ(1));
	t47 = sin(qJ(1));
	t46 = pkin(5) + qJ(2);
	t45 = pkin(7) + qJ(3);
	t44 = cos(t45);
	t43 = sin(t45);
	t42 = cos(pkin(7)) * pkin(2) + pkin(1);
	t1 = [t48 * t44, -t48 * t43, t47, t48 * t42 + t46 * t47 + 0; t47 * t44, -t47 * t43, -t48, t47 * t42 - t48 * t46 + 0; t43, t44, 0, sin(pkin(7)) * pkin(2) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:44:29
	% EndTime: 2020-11-04 19:44:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t54 = sin(qJ(4));
	t55 = sin(qJ(1));
	t62 = t55 * t54;
	t56 = cos(qJ(4));
	t61 = t55 * t56;
	t57 = cos(qJ(1));
	t60 = t57 * t54;
	t59 = t57 * t56;
	t52 = pkin(7) + qJ(3);
	t50 = sin(t52);
	t51 = cos(t52);
	t58 = pkin(3) * t51 + pkin(6) * t50 + cos(pkin(7)) * pkin(2) + pkin(1);
	t53 = pkin(5) + qJ(2);
	t1 = [t51 * t59 + t62, -t51 * t60 + t61, t57 * t50, t53 * t55 + t58 * t57 + 0; t51 * t61 - t60, -t51 * t62 - t59, t55 * t50, -t57 * t53 + t58 * t55 + 0; t50 * t56, -t50 * t54, -t51, t50 * pkin(3) - t51 * pkin(6) + sin(pkin(7)) * pkin(2) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end