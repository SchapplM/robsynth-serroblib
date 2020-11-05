% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:42
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:42:52
	% EndTime: 2020-11-04 20:42:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:42:52
	% EndTime: 2020-11-04 20:42:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t46 = cos(qJ(1));
	t45 = sin(qJ(1));
	t1 = [t46, -t45, 0, 0; t45, t46, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:42:52
	% EndTime: 2020-11-04 20:42:52
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t50 = cos(qJ(1));
	t49 = cos(qJ(2));
	t48 = sin(qJ(1));
	t47 = sin(qJ(2));
	t1 = [t50 * t49, -t50 * t47, t48, t50 * pkin(1) + t48 * pkin(6) + 0; t48 * t49, -t48 * t47, -t50, t48 * pkin(1) - t50 * pkin(6) + 0; t47, t49, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:42:52
	% EndTime: 2020-11-04 20:42:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t57 = pkin(7) + pkin(6);
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t54 = qJ(2) + qJ(3);
	t53 = cos(t54);
	t52 = sin(t54);
	t51 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t56 * t53, -t56 * t52, t55, t56 * t51 + t57 * t55 + 0; t55 * t53, -t55 * t52, -t56, t55 * t51 - t56 * t57 + 0; t52, t53, 0, sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:42:52
	% EndTime: 2020-11-04 20:42:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->15), mult. (14->12), div. (0->0), fcn. (22->8), ass. (0->9)
	t63 = qJ(2) + qJ(3);
	t65 = cos(qJ(1));
	t64 = sin(qJ(1));
	t62 = -qJ(4) - pkin(7) - pkin(6);
	t61 = pkin(9) + t63;
	t60 = cos(t61);
	t59 = sin(t61);
	t58 = pkin(3) * cos(t63) + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t65 * t60, -t65 * t59, t64, t65 * t58 - t64 * t62 + 0; t64 * t60, -t64 * t59, -t65, t64 * t58 + t65 * t62 + 0; t59, t60, 0, pkin(3) * sin(t63) + sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:42:52
	% EndTime: 2020-11-04 20:42:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (60->22), mult. (36->24), div. (0->0), fcn. (49->10), ass. (0->15)
	t72 = sin(qJ(5));
	t73 = sin(qJ(1));
	t80 = t73 * t72;
	t74 = cos(qJ(5));
	t79 = t73 * t74;
	t75 = cos(qJ(1));
	t78 = t75 * t72;
	t77 = t75 * t74;
	t71 = qJ(2) + qJ(3);
	t69 = pkin(9) + t71;
	t67 = sin(t69);
	t68 = cos(t69);
	t76 = pkin(4) * t68 + pkin(8) * t67 + pkin(3) * cos(t71) + cos(qJ(2)) * pkin(2) + pkin(1);
	t70 = -qJ(4) - pkin(7) - pkin(6);
	t1 = [t68 * t77 + t80, -t68 * t78 + t79, t75 * t67, -t73 * t70 + t76 * t75 + 0; t68 * t79 - t78, -t68 * t80 - t77, t73 * t67, t75 * t70 + t76 * t73 + 0; t67 * t74, -t67 * t72, -t68, t67 * pkin(4) - t68 * pkin(8) + pkin(3) * sin(t71) + sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end