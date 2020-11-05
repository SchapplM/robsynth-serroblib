% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:30
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:30:35
	% EndTime: 2020-11-04 20:30:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:30:35
	% EndTime: 2020-11-04 20:30:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t46 = cos(qJ(1));
	t45 = sin(qJ(1));
	t1 = [t46, -t45, 0, 0; t45, t46, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:30:35
	% EndTime: 2020-11-04 20:30:35
	% DurationCPUTime: 0.02s
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
	% StartTime: 2020-11-04 20:30:35
	% EndTime: 2020-11-04 20:30:35
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t55 = -qJ(3) - pkin(6);
	t54 = qJ(2) + pkin(8);
	t53 = cos(t54);
	t52 = sin(t54);
	t51 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t57 * t53, -t57 * t52, t56, t57 * t51 - t55 * t56 + 0; t56 * t53, -t56 * t52, -t57, t56 * t51 + t57 * t55 + 0; t52, t53, 0, sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:30:35
	% EndTime: 2020-11-04 20:30:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->17), mult. (23->17), div. (0->0), fcn. (31->8), ass. (0->11)
	t67 = cos(qJ(1));
	t66 = sin(qJ(1));
	t65 = sin(qJ(2));
	t64 = -qJ(3) - pkin(6);
	t63 = cos(pkin(8));
	t62 = sin(pkin(8));
	t61 = qJ(2) + pkin(8);
	t60 = cos(t61);
	t59 = sin(t61);
	t58 = (pkin(3) * t63 + qJ(4) * t62 + pkin(2)) * cos(qJ(2)) + (-t62 * pkin(3) + qJ(4) * t63) * t65 + pkin(1);
	t1 = [t67 * t60, t66, t67 * t59, t58 * t67 - t64 * t66 + 0; t66 * t60, -t67, t66 * t59, t58 * t66 + t67 * t64 + 0; t59, 0, -t60, t65 * pkin(2) + t59 * pkin(3) - t60 * qJ(4) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:30:35
	% EndTime: 2020-11-04 20:30:35
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (48->21), mult. (35->21), div. (0->0), fcn. (49->10), ass. (0->16)
	t84 = cos(qJ(5));
	t83 = sin(qJ(5));
	t82 = pkin(3) + pkin(4);
	t81 = cos(qJ(1));
	t80 = sin(qJ(1));
	t79 = sin(qJ(2));
	t78 = cos(pkin(8));
	t77 = sin(pkin(8));
	t76 = qJ(2) + pkin(8);
	t75 = qJ(3) + pkin(6) - pkin(7);
	t74 = cos(t76);
	t73 = sin(t76);
	t70 = t73 * t84 - t74 * t83;
	t69 = -t73 * t83 - t74 * t84;
	t68 = (qJ(4) * t77 + t82 * t78 + pkin(2)) * cos(qJ(2)) + (qJ(4) * t78 - t77 * t82) * t79 + pkin(1);
	t1 = [-t81 * t69, t81 * t70, -t80, t68 * t81 + t75 * t80 + 0; -t80 * t69, t80 * t70, t81, t68 * t80 - t75 * t81 + 0; t70, t69, 0, t79 * pkin(2) - t74 * qJ(4) + t82 * t73 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end