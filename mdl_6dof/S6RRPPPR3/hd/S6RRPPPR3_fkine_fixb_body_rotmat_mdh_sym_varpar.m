% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPPR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:59
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:24
	% EndTime: 2020-11-04 21:59:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:24
	% EndTime: 2020-11-04 21:59:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t45 = cos(qJ(1));
	t44 = sin(qJ(1));
	t1 = [t45, -t44, 0, 0; t44, t45, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:24
	% EndTime: 2020-11-04 21:59:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t49 = cos(qJ(1));
	t48 = cos(qJ(2));
	t47 = sin(qJ(1));
	t46 = sin(qJ(2));
	t1 = [t49 * t48, -t49 * t46, t47, t49 * pkin(1) + t47 * pkin(7) + 0; t47 * t48, -t47 * t46, -t49, t47 * pkin(1) - t49 * pkin(7) + 0; t46, t48, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:24
	% EndTime: 2020-11-04 21:59:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->11), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t54 = cos(qJ(1));
	t53 = cos(qJ(2));
	t52 = sin(qJ(1));
	t51 = sin(qJ(2));
	t50 = t53 * pkin(2) + t51 * qJ(3) + pkin(1);
	t1 = [t54 * t53, t52, t54 * t51, t52 * pkin(7) + t50 * t54 + 0; t52 * t53, -t54, t52 * t51, -t54 * pkin(7) + t50 * t52 + 0; t51, 0, -t53, t51 * pkin(2) - t53 * qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:24
	% EndTime: 2020-11-04 21:59:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->16), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->8)
	t61 = pkin(2) + pkin(3);
	t60 = cos(qJ(1));
	t59 = cos(qJ(2));
	t58 = sin(qJ(1));
	t57 = sin(qJ(2));
	t56 = pkin(7) - qJ(4);
	t55 = t57 * qJ(3) + t61 * t59 + pkin(1);
	t1 = [t60 * t57, -t60 * t59, -t58, t55 * t60 + t56 * t58 + 0; t58 * t57, -t58 * t59, t60, t55 * t58 - t56 * t60 + 0; -t59, -t57, 0, -t59 * qJ(3) + t61 * t57 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:24
	% EndTime: 2020-11-04 21:59:24
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (27->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->15)
	t64 = sin(pkin(9));
	t69 = sin(qJ(1));
	t75 = t69 * t64;
	t65 = cos(pkin(9));
	t74 = t69 * t65;
	t71 = cos(qJ(1));
	t73 = t71 * t64;
	t72 = t71 * t65;
	t70 = cos(qJ(2));
	t68 = sin(qJ(2));
	t67 = pkin(4) + qJ(3);
	t66 = pkin(7) - qJ(4);
	t63 = pkin(2) + pkin(3) + qJ(5);
	t62 = t63 * t70 + t67 * t68 + pkin(1);
	t1 = [t68 * t72 - t75, -t68 * t73 - t74, t71 * t70, t62 * t71 + t66 * t69 + 0; t68 * t74 + t73, -t68 * t75 + t72, t69 * t70, t62 * t69 - t66 * t71 + 0; -t70 * t65, t70 * t64, t68, t63 * t68 - t67 * t70 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:24
	% EndTime: 2020-11-04 21:59:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (45->20), mult. (31->22), div. (0->0), fcn. (44->8), ass. (0->16)
	t82 = pkin(9) + qJ(6);
	t79 = sin(t82);
	t84 = sin(qJ(1));
	t90 = t84 * t79;
	t80 = cos(t82);
	t89 = t84 * t80;
	t86 = cos(qJ(1));
	t88 = t86 * t79;
	t87 = t86 * t80;
	t85 = cos(qJ(2));
	t83 = sin(qJ(2));
	t81 = qJ(5) + pkin(2) + pkin(3) + pkin(8);
	t78 = cos(pkin(9)) * pkin(5) + qJ(3) + pkin(4);
	t77 = sin(pkin(9)) * pkin(5) + qJ(4) - pkin(7);
	t76 = t78 * t83 + t81 * t85 + pkin(1);
	t1 = [t83 * t87 - t90, -t83 * t88 - t89, t86 * t85, t76 * t86 - t77 * t84 + 0; t83 * t89 + t88, -t83 * t90 + t87, t84 * t85, t76 * t84 + t77 * t86 + 0; -t85 * t80, t85 * t79, t83, -t78 * t85 + t81 * t83 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end