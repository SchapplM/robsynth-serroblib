% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for the body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:20:30
	% EndTime: 2022-01-23 09:20:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:20:30
	% EndTime: 2022-01-23 09:20:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t45 = cos(qJ(1));
	t44 = sin(qJ(1));
	t1 = [t45, -t44, 0, 0; t44, t45, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:20:30
	% EndTime: 2022-01-23 09:20:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t48 = qJ(1) + pkin(8);
	t47 = cos(t48);
	t46 = sin(t48);
	t1 = [t47, -t46, 0, cos(qJ(1)) * pkin(1) + 0; t46, t47, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:20:30
	% EndTime: 2022-01-23 09:20:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->10), mult. (4->4), div. (0->0), fcn. (8->6), ass. (0->5)
	t52 = qJ(1) + pkin(8);
	t51 = qJ(3) + t52;
	t50 = cos(t51);
	t49 = sin(t51);
	t1 = [t50, -t49, 0, pkin(2) * cos(t52) + cos(qJ(1)) * pkin(1) + 0; t49, t50, 0, pkin(2) * sin(t52) + sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, pkin(6) + qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:20:30
	% EndTime: 2022-01-23 09:20:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (36->16), mult. (12->12), div. (0->0), fcn. (20->8), ass. (0->7)
	t56 = qJ(1) + pkin(8);
	t58 = cos(pkin(9));
	t57 = sin(pkin(9));
	t55 = qJ(3) + t56;
	t54 = cos(t55);
	t53 = sin(t55);
	t1 = [t54 * t58, -t54 * t57, t53, t54 * pkin(3) + t53 * qJ(4) + pkin(2) * cos(t56) + cos(qJ(1)) * pkin(1) + 0; t53 * t58, -t53 * t57, -t54, t53 * pkin(3) - t54 * qJ(4) + pkin(2) * sin(t56) + sin(qJ(1)) * pkin(1) + 0; t57, t58, 0, pkin(6) + qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:20:30
	% EndTime: 2022-01-23 09:20:30
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (61->23), mult. (34->26), div. (0->0), fcn. (47->10), ass. (0->12)
	t64 = cos(pkin(9));
	t65 = sin(qJ(5));
	t69 = t64 * t65;
	t66 = cos(qJ(5));
	t68 = t64 * t66;
	t62 = qJ(1) + pkin(8);
	t63 = sin(pkin(9));
	t67 = pkin(4) * t64 + pkin(7) * t63 + pkin(3);
	t61 = qJ(3) + t62;
	t60 = cos(t61);
	t59 = sin(t61);
	t1 = [t59 * t65 + t60 * t68, t59 * t66 - t60 * t69, t60 * t63, pkin(2) * cos(t62) + cos(qJ(1)) * pkin(1) + t59 * qJ(4) + 0 + t67 * t60; t59 * t68 - t60 * t65, -t59 * t69 - t60 * t66, t59 * t63, pkin(2) * sin(t62) + sin(qJ(1)) * pkin(1) - t60 * qJ(4) + 0 + t67 * t59; t63 * t66, -t63 * t65, -t64, t63 * pkin(4) - t64 * pkin(7) + pkin(5) + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end