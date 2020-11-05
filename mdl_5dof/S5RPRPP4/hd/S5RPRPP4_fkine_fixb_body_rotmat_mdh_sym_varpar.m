% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRPP4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:17
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRPP4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPP4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:43
	% EndTime: 2020-11-04 20:17:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:43
	% EndTime: 2020-11-04 20:17:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t38 = cos(qJ(1));
	t37 = sin(qJ(1));
	t1 = [t38, -t37, 0, 0; t37, t38, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:43
	% EndTime: 2020-11-04 20:17:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t40 = cos(qJ(1));
	t39 = sin(qJ(1));
	t1 = [0, -t40, t39, t40 * pkin(1) + t39 * qJ(2) + 0; 0, -t39, -t40, t39 * pkin(1) - t40 * qJ(2) + 0; 1, 0, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:43
	% EndTime: 2020-11-04 20:17:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t45 = pkin(1) + pkin(6);
	t44 = cos(qJ(1));
	t43 = cos(qJ(3));
	t42 = sin(qJ(1));
	t41 = sin(qJ(3));
	t1 = [t42 * t41, t42 * t43, t44, t42 * qJ(2) + t45 * t44 + 0; -t44 * t41, -t44 * t43, t42, -t44 * qJ(2) + t45 * t42 + 0; t43, -t41, 0, pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:43
	% EndTime: 2020-11-04 20:17:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->13), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t52 = cos(qJ(1));
	t51 = sin(qJ(1));
	t50 = qJ(3) + pkin(7);
	t49 = pkin(1) + pkin(6) + qJ(4);
	t48 = cos(t50);
	t47 = sin(t50);
	t46 = sin(qJ(3)) * pkin(3) + qJ(2);
	t1 = [t51 * t47, t51 * t48, t52, t46 * t51 + t49 * t52 + 0; -t52 * t47, -t52 * t48, t51, -t46 * t52 + t49 * t51 + 0; t48, -t47, 0, cos(qJ(3)) * pkin(3) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:17:43
	% EndTime: 2020-11-04 20:17:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (33->18), mult. (23->17), div. (0->0), fcn. (31->8), ass. (0->11)
	t59 = sin(pkin(7));
	t60 = cos(pkin(7));
	t63 = cos(qJ(3));
	t65 = (pkin(4) * t60 + qJ(5) * t59 + pkin(3)) * sin(qJ(3)) - (-t59 * pkin(4) + qJ(5) * t60) * t63 + qJ(2);
	t64 = cos(qJ(1));
	t62 = sin(qJ(1));
	t58 = qJ(3) + pkin(7);
	t57 = pkin(1) + pkin(6) + qJ(4);
	t56 = cos(t58);
	t55 = sin(t58);
	t1 = [t62 * t55, t64, -t62 * t56, t57 * t64 + t65 * t62 + 0; -t64 * t55, t62, t64 * t56, t57 * t62 - t65 * t64 + 0; t56, 0, t55, t63 * pkin(3) + t56 * pkin(4) + t55 * qJ(5) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end