% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPR9 (for one body)
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
% Datum: 2020-11-04 20:43
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRPR9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:43:39
	% EndTime: 2020-11-04 20:43:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:43:39
	% EndTime: 2020-11-04 20:43:39
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t44 = cos(qJ(1));
	t43 = sin(qJ(1));
	t1 = [t44, -t43, 0, 0; t43, t44, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:43:39
	% EndTime: 2020-11-04 20:43:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t48 = cos(qJ(1));
	t47 = cos(qJ(2));
	t46 = sin(qJ(1));
	t45 = sin(qJ(2));
	t1 = [t48 * t47, -t48 * t45, t46, t48 * pkin(1) + t46 * pkin(6) + 0; t46 * t47, -t46 * t45, -t48, t46 * pkin(1) - t48 * pkin(6) + 0; t45, t47, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:43:39
	% EndTime: 2020-11-04 20:43:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t52 = sin(qJ(1));
	t54 = cos(qJ(2));
	t58 = t52 * t54;
	t50 = sin(qJ(3));
	t55 = cos(qJ(1));
	t57 = t55 * t50;
	t53 = cos(qJ(3));
	t56 = t55 * t53;
	t51 = sin(qJ(2));
	t49 = t54 * pkin(2) + t51 * pkin(7) + pkin(1);
	t1 = [t52 * t50 + t54 * t56, t52 * t53 - t54 * t57, t55 * t51, t52 * pkin(6) + t49 * t55 + 0; t53 * t58 - t57, -t50 * t58 - t56, t52 * t51, -t55 * pkin(6) + t49 * t52 + 0; t51 * t53, -t51 * t50, -t54, t51 * pkin(2) - t54 * pkin(7) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:43:39
	% EndTime: 2020-11-04 20:43:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (31->23), div. (0->0), fcn. (44->8), ass. (0->15)
	t67 = sin(qJ(1));
	t68 = cos(qJ(2));
	t72 = t67 * t68;
	t64 = qJ(3) + pkin(9);
	t62 = sin(t64);
	t69 = cos(qJ(1));
	t71 = t69 * t62;
	t63 = cos(t64);
	t70 = t69 * t63;
	t66 = sin(qJ(2));
	t65 = qJ(4) + pkin(7);
	t61 = cos(qJ(3)) * pkin(3) + pkin(2);
	t60 = sin(qJ(3)) * pkin(3) + pkin(6);
	t59 = t61 * t68 + t65 * t66 + pkin(1);
	t1 = [t67 * t62 + t68 * t70, t67 * t63 - t68 * t71, t69 * t66, t59 * t69 + t60 * t67 + 0; t63 * t72 - t71, -t62 * t72 - t70, t67 * t66, t59 * t67 - t60 * t69 + 0; t66 * t63, -t66 * t62, -t68, t66 * t61 - t68 * t65 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:43:39
	% EndTime: 2020-11-04 20:43:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->25), mult. (44->29), div. (0->0), fcn. (57->11), ass. (0->19)
	t91 = sin(pkin(9)) * pkin(4);
	t82 = sin(qJ(1));
	t84 = cos(qJ(2));
	t90 = t82 * t84;
	t87 = qJ(3) + pkin(9);
	t77 = qJ(5) + t87;
	t74 = sin(t77);
	t85 = cos(qJ(1));
	t89 = t85 * t74;
	t75 = cos(t77);
	t88 = t85 * t75;
	t76 = cos(pkin(9)) * pkin(4) + pkin(3);
	t80 = sin(qJ(3));
	t83 = cos(qJ(3));
	t86 = t76 * t80 + t83 * t91 + pkin(6);
	t81 = sin(qJ(2));
	t78 = pkin(8) + qJ(4) + pkin(7);
	t73 = (t76 * t83 - t80 * t91 + pkin(2)) * t84 + t78 * t81 + pkin(1);
	t1 = [t82 * t74 + t84 * t88, t82 * t75 - t84 * t89, t85 * t81, t73 * t85 + t86 * t82 + 0; t75 * t90 - t89, -t74 * t90 - t88, t82 * t81, t73 * t82 - t86 * t85 + 0; t81 * t75, -t81 * t74, -t84, t81 * (pkin(4) * cos(t87) + t83 * pkin(3) + pkin(2)) - t84 * t78 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end