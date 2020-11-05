% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRR12 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:43
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPRR12_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:43:07
	% EndTime: 2020-11-04 21:43:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:43:07
	% EndTime: 2020-11-04 21:43:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t42 = cos(qJ(1));
	t41 = sin(qJ(1));
	t1 = [t42, -t41, 0, 0; t41, t42, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:43:07
	% EndTime: 2020-11-04 21:43:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t44 = cos(qJ(1));
	t43 = sin(qJ(1));
	t1 = [0, -t44, t43, t44 * pkin(1) + t43 * qJ(2) + 0; 0, -t43, -t44, t43 * pkin(1) - t44 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:43:07
	% EndTime: 2020-11-04 21:43:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t49 = pkin(1) + pkin(7);
	t48 = cos(qJ(1));
	t47 = cos(qJ(3));
	t46 = sin(qJ(1));
	t45 = sin(qJ(3));
	t1 = [t46 * t45, t46 * t47, t48, t46 * qJ(2) + t49 * t48 + 0; -t48 * t45, -t48 * t47, t46, -t48 * qJ(2) + t49 * t46 + 0; t47, -t45, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:43:07
	% EndTime: 2020-11-04 21:43:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (17->14), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->7)
	t55 = pkin(1) + pkin(7);
	t54 = cos(qJ(1));
	t53 = cos(qJ(3));
	t52 = sin(qJ(1));
	t51 = sin(qJ(3));
	t50 = -t51 * pkin(3) + t53 * qJ(4) - qJ(2);
	t1 = [t54, -t52 * t51, -t52 * t53, -t50 * t52 + t55 * t54 + 0; t52, t54 * t51, t54 * t53, t50 * t54 + t55 * t52 + 0; 0, -t53, t51, t53 * pkin(3) + t51 * qJ(4) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:43:07
	% EndTime: 2020-11-04 21:43:07
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (24->17), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->13)
	t59 = sin(qJ(1));
	t61 = cos(qJ(3));
	t67 = t59 * t61;
	t57 = sin(qJ(5));
	t62 = cos(qJ(1));
	t66 = t62 * t57;
	t60 = cos(qJ(5));
	t65 = t62 * t60;
	t58 = sin(qJ(3));
	t63 = pkin(3) + pkin(8);
	t64 = t61 * qJ(4) - t63 * t58 - qJ(2);
	t56 = pkin(1) + pkin(4) + pkin(7);
	t1 = [-t57 * t67 + t65, -t60 * t67 - t66, t59 * t58, t56 * t62 - t64 * t59 + 0; t59 * t60 + t61 * t66, -t59 * t57 + t61 * t65, -t62 * t58, t56 * t59 + t64 * t62 + 0; t58 * t57, t58 * t60, t61, t58 * qJ(4) + t63 * t61 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:43:07
	% EndTime: 2020-11-04 21:43:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (42->21), mult. (31->23), div. (0->0), fcn. (44->8), ass. (0->15)
	t75 = sin(qJ(1));
	t76 = cos(qJ(3));
	t81 = t75 * t76;
	t73 = qJ(5) + qJ(6);
	t70 = sin(t73);
	t77 = cos(qJ(1));
	t80 = t77 * t70;
	t71 = cos(t73);
	t79 = t77 * t71;
	t69 = sin(qJ(5)) * pkin(5) + qJ(4);
	t72 = pkin(3) + pkin(8) + pkin(9);
	t74 = sin(qJ(3));
	t78 = t69 * t76 - t72 * t74 - qJ(2);
	t68 = cos(qJ(5)) * pkin(5) + pkin(1) + pkin(4) + pkin(7);
	t1 = [-t70 * t81 + t79, -t71 * t81 - t80, t75 * t74, t68 * t77 - t78 * t75 + 0; t75 * t71 + t76 * t80, -t75 * t70 + t76 * t79, -t77 * t74, t68 * t75 + t78 * t77 + 0; t74 * t70, t74 * t71, t76, t69 * t74 + t72 * t76 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end