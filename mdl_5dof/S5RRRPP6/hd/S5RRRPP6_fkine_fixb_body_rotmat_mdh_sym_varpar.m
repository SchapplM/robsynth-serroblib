% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPP6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:41
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRPP6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:41:01
	% EndTime: 2020-11-04 20:41:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:41:01
	% EndTime: 2020-11-04 20:41:01
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t49 = cos(qJ(1));
	t48 = sin(qJ(1));
	t1 = [t49, -t48, 0, 0; t48, t49, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:41:01
	% EndTime: 2020-11-04 20:41:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t53 = cos(qJ(1));
	t52 = cos(qJ(2));
	t51 = sin(qJ(1));
	t50 = sin(qJ(2));
	t1 = [t53 * t52, -t53 * t50, t51, t53 * pkin(1) + t51 * pkin(6) + 0; t51 * t52, -t51 * t50, -t53, t51 * pkin(1) - t53 * pkin(6) + 0; t50, t52, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:41:01
	% EndTime: 2020-11-04 20:41:01
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t57 = sin(qJ(1));
	t59 = cos(qJ(2));
	t63 = t57 * t59;
	t55 = sin(qJ(3));
	t60 = cos(qJ(1));
	t62 = t60 * t55;
	t58 = cos(qJ(3));
	t61 = t60 * t58;
	t56 = sin(qJ(2));
	t54 = t59 * pkin(2) + t56 * pkin(7) + pkin(1);
	t1 = [t57 * t55 + t59 * t61, t57 * t58 - t59 * t62, t60 * t56, t57 * pkin(6) + t54 * t60 + 0; t58 * t63 - t62, -t55 * t63 - t61, t57 * t56, -t60 * pkin(6) + t54 * t57 + 0; t56 * t58, -t56 * t55, -t59, t56 * pkin(2) - t59 * pkin(7) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:41:01
	% EndTime: 2020-11-04 20:41:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (31->23), div. (0->0), fcn. (44->8), ass. (0->15)
	t72 = sin(qJ(1));
	t73 = cos(qJ(2));
	t77 = t72 * t73;
	t69 = qJ(3) + pkin(8);
	t67 = sin(t69);
	t74 = cos(qJ(1));
	t76 = t74 * t67;
	t68 = cos(t69);
	t75 = t74 * t68;
	t71 = sin(qJ(2));
	t70 = qJ(4) + pkin(7);
	t66 = cos(qJ(3)) * pkin(3) + pkin(2);
	t65 = sin(qJ(3)) * pkin(3) + pkin(6);
	t64 = t66 * t73 + t70 * t71 + pkin(1);
	t1 = [t72 * t67 + t73 * t75, t72 * t68 - t73 * t76, t74 * t71, t64 * t74 + t65 * t72 + 0; t68 * t77 - t76, -t67 * t77 - t75, t72 * t71, t64 * t72 - t65 * t74 + 0; t71 * t68, -t71 * t67, -t73, t71 * t66 - t73 * t70 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:41:02
	% EndTime: 2020-11-04 20:41:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (54->23), mult. (56->29), div. (0->0), fcn. (69->10), ass. (0->21)
	t90 = sin(qJ(1));
	t92 = cos(qJ(2));
	t97 = t90 * t92;
	t84 = qJ(3) + pkin(8);
	t82 = sin(t84);
	t93 = cos(qJ(1));
	t96 = t93 * t82;
	t83 = cos(t84);
	t95 = t93 * t83;
	t85 = sin(pkin(8));
	t86 = cos(pkin(8));
	t80 = pkin(4) * t86 + qJ(5) * t85 + pkin(3);
	t81 = -t85 * pkin(4) + qJ(5) * t86;
	t88 = sin(qJ(3));
	t91 = cos(qJ(3));
	t94 = t80 * t88 - t81 * t91 + pkin(6);
	t89 = sin(qJ(2));
	t87 = qJ(4) + pkin(7);
	t79 = t80 * t91 + t81 * t88 + pkin(2);
	t78 = t79 * t92 + t87 * t89 + pkin(1);
	t1 = [t90 * t82 + t92 * t95, t93 * t89, -t90 * t83 + t92 * t96, t78 * t93 + t94 * t90 + 0; t83 * t97 - t96, t90 * t89, t82 * t97 + t95, t78 * t90 - t94 * t93 + 0; t89 * t83, -t92, t89 * t82, t79 * t89 - t92 * t87 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end