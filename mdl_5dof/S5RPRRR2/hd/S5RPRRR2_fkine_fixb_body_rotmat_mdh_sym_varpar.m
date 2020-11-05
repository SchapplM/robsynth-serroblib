% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:26
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:26:09
	% EndTime: 2020-11-04 20:26:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:26:09
	% EndTime: 2020-11-04 20:26:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t38 = cos(qJ(1));
	t37 = sin(qJ(1));
	t1 = [t38, -t37, 0, 0; t37, t38, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:26:09
	% EndTime: 2020-11-04 20:26:09
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t42 = cos(qJ(1));
	t41 = sin(qJ(1));
	t40 = cos(pkin(9));
	t39 = sin(pkin(9));
	t1 = [t42 * t40, -t42 * t39, t41, t42 * pkin(1) + t41 * qJ(2) + 0; t41 * t40, -t41 * t39, -t42, t41 * pkin(1) - t42 * qJ(2) + 0; t39, t40, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:26:09
	% EndTime: 2020-11-04 20:26:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t49 = cos(qJ(1));
	t48 = sin(qJ(1));
	t47 = pkin(6) + qJ(2);
	t46 = pkin(9) + qJ(3);
	t45 = cos(t46);
	t44 = sin(t46);
	t43 = cos(pkin(9)) * pkin(2) + pkin(1);
	t1 = [t49 * t45, -t49 * t44, t48, t49 * t43 + t47 * t48 + 0; t48 * t45, -t48 * t44, -t49, t48 * t43 - t49 * t47 + 0; t44, t45, 0, sin(pkin(9)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:26:09
	% EndTime: 2020-11-04 20:26:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->15), mult. (14->12), div. (0->0), fcn. (22->8), ass. (0->9)
	t55 = pkin(9) + qJ(3);
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t54 = -pkin(7) - pkin(6) - qJ(2);
	t53 = qJ(4) + t55;
	t52 = cos(t53);
	t51 = sin(t53);
	t50 = pkin(3) * cos(t55) + cos(pkin(9)) * pkin(2) + pkin(1);
	t1 = [t57 * t52, -t57 * t51, t56, t57 * t50 - t56 * t54 + 0; t56 * t52, -t56 * t51, -t57, t56 * t50 + t57 * t54 + 0; t51, t52, 0, pkin(3) * sin(t55) + sin(pkin(9)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:26:09
	% EndTime: 2020-11-04 20:26:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (50->18), mult. (17->14), div. (0->0), fcn. (25->10), ass. (0->10)
	t64 = pkin(9) + qJ(3);
	t62 = qJ(4) + t64;
	t66 = cos(qJ(1));
	t65 = sin(qJ(1));
	t63 = -pkin(8) - pkin(7) - pkin(6) - qJ(2);
	t61 = qJ(5) + t62;
	t60 = cos(t61);
	t59 = sin(t61);
	t58 = pkin(4) * cos(t62) + pkin(3) * cos(t64) + cos(pkin(9)) * pkin(2) + pkin(1);
	t1 = [t66 * t60, -t66 * t59, t65, t66 * t58 - t65 * t63 + 0; t65 * t60, -t65 * t59, -t66, t65 * t58 + t66 * t63 + 0; t59, t60, 0, pkin(4) * sin(t62) + pkin(3) * sin(t64) + sin(pkin(9)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end