% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRR13 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:28
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRRR13_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR13_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:28:36
	% EndTime: 2020-11-04 20:28:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:28:36
	% EndTime: 2020-11-04 20:28:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t38 = cos(qJ(1));
	t37 = sin(qJ(1));
	t1 = [t38, -t37, 0, 0; t37, t38, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:28:36
	% EndTime: 2020-11-04 20:28:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t40 = cos(qJ(1));
	t39 = sin(qJ(1));
	t1 = [0, -t40, t39, t40 * pkin(1) + t39 * qJ(2) + 0; 0, -t39, -t40, t39 * pkin(1) - t40 * qJ(2) + 0; 1, 0, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:28:36
	% EndTime: 2020-11-04 20:28:37
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
	% StartTime: 2020-11-04 20:28:37
	% EndTime: 2020-11-04 20:28:37
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (20->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->13)
	t47 = sin(qJ(4));
	t49 = sin(qJ(1));
	t57 = t49 * t47;
	t50 = cos(qJ(4));
	t56 = t49 * t50;
	t52 = cos(qJ(1));
	t55 = t52 * t47;
	t54 = t52 * t50;
	t53 = pkin(1) + pkin(6);
	t51 = cos(qJ(3));
	t48 = sin(qJ(3));
	t46 = -t48 * pkin(3) + t51 * pkin(7) - qJ(2);
	t1 = [t48 * t56 + t55, -t48 * t57 + t54, -t49 * t51, -t46 * t49 + t53 * t52 + 0; -t48 * t54 + t57, t48 * t55 + t56, t52 * t51, t46 * t52 + t53 * t49 + 0; t51 * t50, -t51 * t47, t48, t51 * pkin(3) + t48 * pkin(7) + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:28:37
	% EndTime: 2020-11-04 20:28:37
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (38->21), mult. (31->22), div. (0->0), fcn. (44->8), ass. (0->16)
	t62 = qJ(4) + qJ(5);
	t60 = sin(t62);
	t64 = sin(qJ(1));
	t72 = t64 * t60;
	t61 = cos(t62);
	t71 = t64 * t61;
	t66 = cos(qJ(1));
	t70 = t66 * t60;
	t69 = t66 * t61;
	t59 = cos(qJ(4)) * pkin(4) + pkin(3);
	t63 = sin(qJ(3));
	t65 = cos(qJ(3));
	t67 = pkin(8) + pkin(7);
	t68 = t59 * t63 - t67 * t65 + qJ(2);
	t58 = sin(qJ(4)) * pkin(4) + pkin(1) + pkin(6);
	t1 = [t63 * t71 + t70, -t63 * t72 + t69, -t64 * t65, t58 * t66 + t68 * t64 + 0; -t63 * t69 + t72, t63 * t70 + t71, t66 * t65, t58 * t64 - t68 * t66 + 0; t65 * t61, -t65 * t60, t63, t65 * t59 + t67 * t63 + pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end