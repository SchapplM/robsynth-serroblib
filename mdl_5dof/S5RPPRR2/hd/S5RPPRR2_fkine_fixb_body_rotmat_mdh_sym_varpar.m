% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:13
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:57
	% EndTime: 2020-11-04 20:13:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:57
	% EndTime: 2020-11-04 20:13:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t37 = cos(qJ(1));
	t36 = sin(qJ(1));
	t1 = [t37, -t36, 0, 0; t36, t37, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:57
	% EndTime: 2020-11-04 20:13:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t39 = cos(qJ(1));
	t38 = sin(qJ(1));
	t1 = [0, -t39, t38, t39 * pkin(1) + t38 * qJ(2) + 0; 0, -t38, -t39, t38 * pkin(1) - t39 * qJ(2) + 0; 1, 0, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:57
	% EndTime: 2020-11-04 20:13:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t44 = cos(qJ(1));
	t43 = sin(qJ(1));
	t42 = pkin(1) + qJ(3);
	t41 = cos(pkin(8));
	t40 = sin(pkin(8));
	t1 = [t43 * t40, t43 * t41, t44, t43 * qJ(2) + t42 * t44 + 0; -t44 * t40, -t44 * t41, t43, -t44 * qJ(2) + t42 * t43 + 0; t41, -t40, 0, pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:57
	% EndTime: 2020-11-04 20:13:57
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (22->13), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t51 = cos(qJ(1));
	t50 = sin(qJ(1));
	t49 = pkin(8) + qJ(4);
	t48 = pkin(1) + pkin(6) + qJ(3);
	t47 = cos(t49);
	t46 = sin(t49);
	t45 = sin(pkin(8)) * pkin(3) + qJ(2);
	t1 = [t50 * t46, t50 * t47, t51, t45 * t50 + t48 * t51 + 0; -t51 * t46, -t51 * t47, t50, -t45 * t51 + t48 * t50 + 0; t47, -t46, 0, cos(pkin(8)) * pkin(3) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:13:57
	% EndTime: 2020-11-04 20:13:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (36->16), mult. (16->12), div. (0->0), fcn. (24->8), ass. (0->9)
	t58 = pkin(8) + qJ(4);
	t61 = pkin(4) * sin(t58) + sin(pkin(8)) * pkin(3) + qJ(2);
	t60 = cos(qJ(1));
	t59 = sin(qJ(1));
	t57 = qJ(3) + pkin(1) + pkin(6) + pkin(7);
	t56 = qJ(5) + t58;
	t54 = cos(t56);
	t53 = sin(t56);
	t1 = [t59 * t53, t59 * t54, t60, t57 * t60 + t61 * t59 + 0; -t60 * t53, -t60 * t54, t59, t57 * t59 - t61 * t60 + 0; t54, -t53, 0, pkin(4) * cos(t58) + cos(pkin(8)) * pkin(3) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end