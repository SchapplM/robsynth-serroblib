% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR9 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:50
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:50:50
	% EndTime: 2020-11-04 20:50:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:50:50
	% EndTime: 2020-11-04 20:50:50
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t46 = cos(qJ(1));
	t45 = sin(qJ(1));
	t1 = [t46, -t45, 0, 0; t45, t46, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:50:50
	% EndTime: 2020-11-04 20:50:50
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
	% StartTime: 2020-11-04 20:50:50
	% EndTime: 2020-11-04 20:50:50
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t54 = sin(qJ(1));
	t56 = cos(qJ(2));
	t60 = t54 * t56;
	t52 = sin(qJ(3));
	t57 = cos(qJ(1));
	t59 = t57 * t52;
	t55 = cos(qJ(3));
	t58 = t57 * t55;
	t53 = sin(qJ(2));
	t51 = t56 * pkin(2) + t53 * pkin(7) + pkin(1);
	t1 = [t54 * t52 + t56 * t58, t54 * t55 - t56 * t59, t57 * t53, t54 * pkin(6) + t51 * t57 + 0; t55 * t60 - t59, -t52 * t60 - t58, t54 * t53, -t57 * pkin(6) + t51 * t54 + 0; t53 * t55, -t53 * t52, -t56, t53 * pkin(2) - t56 * pkin(7) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:50:50
	% EndTime: 2020-11-04 20:50:50
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (35->19), mult. (31->23), div. (0->0), fcn. (44->8), ass. (0->15)
	t68 = sin(qJ(1));
	t69 = cos(qJ(2));
	t74 = t68 * t69;
	t66 = qJ(3) + qJ(4);
	t64 = sin(t66);
	t70 = cos(qJ(1));
	t73 = t70 * t64;
	t65 = cos(t66);
	t72 = t70 * t65;
	t71 = pkin(8) + pkin(7);
	t67 = sin(qJ(2));
	t63 = cos(qJ(3)) * pkin(3) + pkin(2);
	t62 = sin(qJ(3)) * pkin(3) + pkin(6);
	t61 = t63 * t69 + t71 * t67 + pkin(1);
	t1 = [t68 * t64 + t69 * t72, t68 * t65 - t69 * t73, t70 * t67, t61 * t70 + t62 * t68 + 0; t65 * t74 - t73, -t64 * t74 - t72, t68 * t67, t61 * t68 - t62 * t70 + 0; t67 * t65, -t67 * t64, -t69, t67 * t63 - t69 * t71 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:50:50
	% EndTime: 2020-11-04 20:50:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->22), mult. (42->25), div. (0->0), fcn. (55->10), ass. (0->16)
	t81 = qJ(3) + qJ(4);
	t90 = pkin(6) + pkin(4) * sin(t81) + sin(qJ(3)) * pkin(3);
	t83 = sin(qJ(1));
	t84 = cos(qJ(2));
	t89 = t83 * t84;
	t79 = qJ(5) + t81;
	t77 = sin(t79);
	t85 = cos(qJ(1));
	t88 = t85 * t77;
	t78 = cos(t79);
	t87 = t85 * t78;
	t75 = pkin(4) * cos(t81) + cos(qJ(3)) * pkin(3) + pkin(2);
	t80 = -pkin(9) - pkin(8) - pkin(7);
	t82 = sin(qJ(2));
	t86 = t75 * t84 - t80 * t82 + pkin(1);
	t1 = [t83 * t77 + t84 * t87, t83 * t78 - t84 * t88, t85 * t82, t90 * t83 + t86 * t85 + 0; t78 * t89 - t88, -t77 * t89 - t87, t83 * t82, t86 * t83 - t90 * t85 + 0; t82 * t78, -t82 * t77, -t84, t82 * t75 + t84 * t80 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end