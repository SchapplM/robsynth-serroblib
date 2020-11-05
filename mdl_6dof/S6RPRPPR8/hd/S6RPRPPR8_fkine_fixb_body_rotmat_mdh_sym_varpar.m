% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPPR8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:35
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPPR8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:35:27
	% EndTime: 2020-11-04 21:35:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:35:27
	% EndTime: 2020-11-04 21:35:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t40 = cos(qJ(1));
	t39 = sin(qJ(1));
	t1 = [t40, -t39, 0, 0; t39, t40, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:35:27
	% EndTime: 2020-11-04 21:35:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t42 = cos(qJ(1));
	t41 = sin(qJ(1));
	t1 = [0, -t42, t41, t42 * pkin(1) + t41 * qJ(2) + 0; 0, -t41, -t42, t41 * pkin(1) - t42 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:35:27
	% EndTime: 2020-11-04 21:35:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t47 = pkin(1) + pkin(7);
	t46 = cos(qJ(1));
	t45 = cos(qJ(3));
	t44 = sin(qJ(1));
	t43 = sin(qJ(3));
	t1 = [t44 * t43, t44 * t45, t46, t44 * qJ(2) + t47 * t46 + 0; -t46 * t43, -t46 * t45, t44, -t46 * qJ(2) + t47 * t44 + 0; t45, -t43, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:35:27
	% EndTime: 2020-11-04 21:35:28
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (16->13), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->7)
	t53 = pkin(1) + pkin(7);
	t52 = cos(qJ(1));
	t51 = cos(qJ(3));
	t50 = sin(qJ(1));
	t49 = sin(qJ(3));
	t48 = -t49 * pkin(3) + t51 * qJ(4) - qJ(2);
	t1 = [t50 * t49, t52, -t50 * t51, -t48 * t50 + t53 * t52 + 0; -t52 * t49, t50, t52 * t51, t48 * t52 + t53 * t50 + 0; t51, 0, t49, t51 * pkin(3) + t49 * qJ(4) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:35:28
	% EndTime: 2020-11-04 21:35:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (24->17), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->8)
	t55 = sin(qJ(3));
	t57 = cos(qJ(3));
	t59 = pkin(3) + pkin(4);
	t60 = t57 * qJ(4) - t59 * t55 - qJ(2);
	t58 = cos(qJ(1));
	t56 = sin(qJ(1));
	t54 = -qJ(5) + pkin(1) + pkin(7);
	t1 = [-t56 * t57, -t56 * t55, -t58, t54 * t58 - t60 * t56 + 0; t58 * t57, t58 * t55, -t56, t54 * t56 + t60 * t58 + 0; t55, -t57, 0, t55 * qJ(4) + t59 * t57 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:35:28
	% EndTime: 2020-11-04 21:35:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->19), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->14)
	t66 = sin(qJ(1));
	t68 = cos(qJ(3));
	t73 = t66 * t68;
	t64 = sin(qJ(6));
	t69 = cos(qJ(1));
	t72 = t69 * t64;
	t67 = cos(qJ(6));
	t71 = t69 * t67;
	t62 = pkin(3) + pkin(4) + pkin(8);
	t63 = -qJ(4) - pkin(5);
	t65 = sin(qJ(3));
	t70 = t62 * t65 + t63 * t68 + qJ(2);
	t61 = -qJ(5) + pkin(1) + pkin(7);
	t1 = [-t67 * t73 - t72, t64 * t73 - t71, t66 * t65, t61 * t69 + t70 * t66 + 0; -t66 * t64 + t68 * t71, -t66 * t67 - t68 * t72, -t69 * t65, t61 * t66 - t70 * t69 + 0; t65 * t67, -t65 * t64, t68, t62 * t68 - t63 * t65 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end