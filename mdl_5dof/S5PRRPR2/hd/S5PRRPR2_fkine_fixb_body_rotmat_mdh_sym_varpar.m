% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRPR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:03
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:03:19
	% EndTime: 2020-11-04 20:03:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:03:19
	% EndTime: 2020-11-04 20:03:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t43 = cos(pkin(8));
	t42 = sin(pkin(8));
	t1 = [t43, -t42, 0, 0; t42, t43, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:03:19
	% EndTime: 2020-11-04 20:03:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t46 = pkin(8) + qJ(2);
	t45 = cos(t46);
	t44 = sin(t46);
	t1 = [t45, -t44, 0, cos(pkin(8)) * pkin(1) + 0; t44, t45, 0, sin(pkin(8)) * pkin(1) + 0; 0, 0, 1, pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:03:19
	% EndTime: 2020-11-04 20:03:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->10), mult. (4->4), div. (0->0), fcn. (8->6), ass. (0->5)
	t50 = pkin(8) + qJ(2);
	t49 = qJ(3) + t50;
	t48 = cos(t49);
	t47 = sin(t49);
	t1 = [t48, -t47, 0, pkin(2) * cos(t50) + cos(pkin(8)) * pkin(1) + 0; t47, t48, 0, pkin(2) * sin(t50) + sin(pkin(8)) * pkin(1) + 0; 0, 0, 1, pkin(6) + pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:03:19
	% EndTime: 2020-11-04 20:03:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (36->16), mult. (12->12), div. (0->0), fcn. (20->8), ass. (0->7)
	t54 = pkin(8) + qJ(2);
	t56 = cos(pkin(9));
	t55 = sin(pkin(9));
	t53 = qJ(3) + t54;
	t52 = cos(t53);
	t51 = sin(t53);
	t1 = [t52 * t56, -t52 * t55, t51, t52 * pkin(3) + t51 * qJ(4) + pkin(2) * cos(t54) + cos(pkin(8)) * pkin(1) + 0; t51 * t56, -t51 * t55, -t52, t51 * pkin(3) - t52 * qJ(4) + pkin(2) * sin(t54) + sin(pkin(8)) * pkin(1) + 0; t55, t56, 0, pkin(6) + pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:03:19
	% EndTime: 2020-11-04 20:03:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (53->23), mult. (30->26), div. (0->0), fcn. (43->10), ass. (0->12)
	t63 = cos(pkin(9));
	t64 = sin(qJ(5));
	t67 = t63 * t64;
	t65 = cos(qJ(5));
	t66 = t63 * t65;
	t61 = pkin(8) + qJ(2);
	t62 = sin(pkin(9));
	t60 = qJ(3) + t61;
	t59 = cos(t60);
	t58 = sin(t60);
	t57 = t63 * pkin(4) + t62 * pkin(7) + pkin(3);
	t1 = [t58 * t64 + t59 * t66, t58 * t65 - t59 * t67, t59 * t62, t57 * t59 + t58 * qJ(4) + cos(pkin(8)) * pkin(1) + pkin(2) * cos(t61) + 0; t58 * t66 - t59 * t64, -t58 * t67 - t59 * t65, t58 * t62, t57 * t58 - t59 * qJ(4) + sin(pkin(8)) * pkin(1) + pkin(2) * sin(t61) + 0; t62 * t65, -t62 * t64, -t63, t62 * pkin(4) - t63 * pkin(7) + pkin(5) + pkin(6) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end