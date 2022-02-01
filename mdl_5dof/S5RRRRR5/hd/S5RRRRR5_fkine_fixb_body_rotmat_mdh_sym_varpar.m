% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for the body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:01:34
	% EndTime: 2022-01-20 12:01:34
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:01:34
	% EndTime: 2022-01-20 12:01:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t38 = cos(qJ(1));
	t37 = sin(qJ(1));
	t1 = [t38, -t37, 0, 0; t37, t38, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:01:34
	% EndTime: 2022-01-20 12:01:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t41 = qJ(1) + qJ(2);
	t40 = cos(t41);
	t39 = sin(t41);
	t1 = [t40, -t39, 0, cos(qJ(1)) * pkin(1) + 0; t39, t40, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:01:34
	% EndTime: 2022-01-20 12:01:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->10), mult. (4->4), div. (0->0), fcn. (8->6), ass. (0->5)
	t45 = qJ(1) + qJ(2);
	t44 = qJ(3) + t45;
	t43 = cos(t44);
	t42 = sin(t44);
	t1 = [t43, -t42, 0, pkin(2) * cos(t45) + cos(qJ(1)) * pkin(1) + 0; t42, t43, 0, pkin(2) * sin(t45) + sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, pkin(7) + pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:01:34
	% EndTime: 2022-01-20 12:01:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (36->16), mult. (12->12), div. (0->0), fcn. (20->8), ass. (0->7)
	t49 = qJ(1) + qJ(2);
	t51 = cos(qJ(4));
	t50 = sin(qJ(4));
	t48 = qJ(3) + t49;
	t47 = cos(t48);
	t46 = sin(t48);
	t1 = [t47 * t51, -t47 * t50, t46, t47 * pkin(3) + t46 * pkin(8) + pkin(2) * cos(t49) + cos(qJ(1)) * pkin(1) + 0; t46 * t51, -t46 * t50, -t47, t46 * pkin(3) - t47 * pkin(8) + pkin(2) * sin(t49) + sin(qJ(1)) * pkin(1) + 0; t50, t51, 0, pkin(7) + pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:01:34
	% EndTime: 2022-01-20 12:01:34
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (47->20), mult. (15->14), div. (0->0), fcn. (23->10), ass. (0->10)
	t59 = qJ(1) + qJ(2);
	t60 = -pkin(9) - pkin(8);
	t58 = qJ(4) + qJ(5);
	t57 = qJ(3) + t59;
	t56 = cos(t58);
	t55 = sin(t58);
	t54 = cos(qJ(4)) * pkin(4) + pkin(3);
	t53 = cos(t57);
	t52 = sin(t57);
	t1 = [t53 * t56, -t53 * t55, t52, t53 * t54 - t52 * t60 + pkin(2) * cos(t59) + cos(qJ(1)) * pkin(1) + 0; t52 * t56, -t52 * t55, -t53, t52 * t54 + t53 * t60 + pkin(2) * sin(t59) + sin(qJ(1)) * pkin(1) + 0; t55, t56, 0, sin(qJ(4)) * pkin(4) + pkin(7) + pkin(6) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end