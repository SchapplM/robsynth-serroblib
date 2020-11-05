% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPPRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:24
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPPPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:24:19
	% EndTime: 2020-11-04 21:24:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:24:19
	% EndTime: 2020-11-04 21:24:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t44 = cos(qJ(1));
	t43 = sin(qJ(1));
	t1 = [t44, -t43, 0, 0; t43, t44, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:24:19
	% EndTime: 2020-11-04 21:24:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t46 = cos(qJ(1));
	t45 = sin(qJ(1));
	t1 = [0, -t46, t45, t46 * pkin(1) + t45 * qJ(2) + 0; 0, -t45, -t46, t45 * pkin(1) - t46 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:24:19
	% EndTime: 2020-11-04 21:24:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->4)
	t49 = cos(qJ(1));
	t48 = sin(qJ(1));
	t47 = pkin(1) + qJ(3);
	t1 = [t48, 0, t49, t48 * qJ(2) + t47 * t49 + 0; -t49, 0, t48, -t49 * qJ(2) + t47 * t48 + 0; 0, -1, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:24:19
	% EndTime: 2020-11-04 21:24:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (15->12), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->9)
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t55 = pkin(1) + qJ(3);
	t54 = pkin(3) + qJ(2);
	t53 = cos(pkin(9));
	t52 = sin(pkin(9));
	t51 = t57 * t52 + t56 * t53;
	t50 = -t56 * t52 + t57 * t53;
	t1 = [t51, t50, 0, t54 * t56 + t55 * t57 + 0; -t50, t51, 0, -t54 * t57 + t55 * t56 + 0; 0, 0, 1, qJ(4) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:24:19
	% EndTime: 2020-11-04 21:24:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (26->19), mult. (28->16), div. (0->0), fcn. (42->6), ass. (0->11)
	t68 = cos(qJ(1));
	t67 = cos(qJ(5));
	t66 = sin(qJ(1));
	t65 = sin(qJ(5));
	t64 = cos(pkin(9));
	t63 = sin(pkin(9));
	t61 = t64 * pkin(4) + t63 * pkin(7) + pkin(3) + qJ(2);
	t60 = t63 * pkin(4) - t64 * pkin(7) + pkin(1) + qJ(3);
	t59 = t68 * t63 + t66 * t64;
	t58 = -t66 * t63 + t68 * t64;
	t1 = [t59 * t67, -t59 * t65, -t58, t60 * t68 + t61 * t66 + 0; -t58 * t67, t58 * t65, -t59, t60 * t66 - t61 * t68 + 0; t65, t67, 0, qJ(4) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:24:19
	% EndTime: 2020-11-04 21:24:19
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (46->26), mult. (66->30), div. (0->0), fcn. (89->8), ass. (0->16)
	t84 = sin(qJ(1));
	t76 = sin(qJ(6));
	t79 = cos(qJ(5));
	t83 = t76 * t79;
	t78 = cos(qJ(6));
	t82 = t78 * t79;
	t77 = sin(qJ(5));
	t81 = pkin(5) * t79 + pkin(8) * t77 + pkin(4);
	t80 = cos(qJ(1));
	t75 = cos(pkin(9));
	t74 = sin(pkin(9));
	t72 = t74 * t80 + t75 * t84;
	t71 = -t74 * t84 + t75 * t80;
	t70 = t74 * pkin(7) + t75 * t81 + pkin(3) + qJ(2);
	t69 = -t75 * pkin(7) + t74 * t81 + pkin(1) + qJ(3);
	t1 = [-t71 * t76 + t72 * t82, -t71 * t78 - t72 * t83, t72 * t77, t69 * t80 + t70 * t84 + 0; -t71 * t82 - t72 * t76, t71 * t83 - t72 * t78, -t71 * t77, t69 * t84 - t70 * t80 + 0; t77 * t78, -t77 * t76, -t79, pkin(5) * t77 - pkin(8) * t79 + pkin(2) + pkin(6) + qJ(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end