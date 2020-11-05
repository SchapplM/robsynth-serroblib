% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRP9 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:46
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRRP9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:58
	% EndTime: 2020-11-04 20:46:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:58
	% EndTime: 2020-11-04 20:46:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t48 = cos(qJ(1));
	t47 = sin(qJ(1));
	t1 = [t48, -t47, 0, 0; t47, t48, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:58
	% EndTime: 2020-11-04 20:46:58
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t52 = cos(qJ(1));
	t51 = cos(qJ(2));
	t50 = sin(qJ(1));
	t49 = sin(qJ(2));
	t1 = [t52 * t51, -t52 * t49, t50, pkin(1) * t52 + pkin(6) * t50 + 0; t50 * t51, -t50 * t49, -t52, pkin(1) * t50 - pkin(6) * t52 + 0; t49, t51, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:58
	% EndTime: 2020-11-04 20:46:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t56 = sin(qJ(1));
	t58 = cos(qJ(2));
	t62 = t56 * t58;
	t54 = sin(qJ(3));
	t59 = cos(qJ(1));
	t61 = t59 * t54;
	t57 = cos(qJ(3));
	t60 = t59 * t57;
	t55 = sin(qJ(2));
	t53 = t58 * pkin(2) + t55 * pkin(7) + pkin(1);
	t1 = [t56 * t54 + t58 * t60, t56 * t57 - t58 * t61, t59 * t55, t56 * pkin(6) + t53 * t59 + 0; t57 * t62 - t61, -t54 * t62 - t60, t56 * t55, -t59 * pkin(6) + t53 * t56 + 0; t55 * t57, -t55 * t54, -t58, t55 * pkin(2) - t58 * pkin(7) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:58
	% EndTime: 2020-11-04 20:46:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (31->23), div. (0->0), fcn. (44->8), ass. (0->15)
	t70 = sin(qJ(1));
	t71 = cos(qJ(2));
	t76 = t70 * t71;
	t68 = qJ(3) + qJ(4);
	t66 = sin(t68);
	t72 = cos(qJ(1));
	t75 = t72 * t66;
	t67 = cos(t68);
	t74 = t72 * t67;
	t73 = pkin(8) + pkin(7);
	t69 = sin(qJ(2));
	t65 = cos(qJ(3)) * pkin(3) + pkin(2);
	t64 = sin(qJ(3)) * pkin(3) + pkin(6);
	t63 = t65 * t71 + t73 * t69 + pkin(1);
	t1 = [t70 * t66 + t71 * t74, t70 * t67 - t71 * t75, t72 * t69, t63 * t72 + t64 * t70 + 0; t67 * t76 - t75, -t66 * t76 - t74, t70 * t69, t63 * t70 - t64 * t72 + 0; t69 * t67, -t69 * t66, -t71, t69 * t65 - t71 * t73 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:58
	% EndTime: 2020-11-04 20:46:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (54->24), mult. (51->29), div. (0->0), fcn. (68->8), ass. (0->19)
	t88 = sin(qJ(1));
	t89 = cos(qJ(2));
	t94 = t88 * t89;
	t86 = qJ(3) + qJ(4);
	t84 = sin(t86);
	t90 = cos(qJ(1));
	t93 = t90 * t84;
	t85 = cos(t86);
	t92 = t90 * t85;
	t91 = pkin(8) + pkin(7);
	t87 = sin(qJ(2));
	t83 = cos(qJ(3)) * pkin(3) + pkin(2);
	t82 = sin(qJ(3)) * pkin(3) + pkin(6);
	t81 = t83 * t89 + t91 * t87 + pkin(1);
	t80 = t88 * t84 + t89 * t92;
	t79 = -t88 * t85 + t89 * t93;
	t78 = t85 * t94 - t93;
	t77 = t84 * t94 + t92;
	t1 = [t80, t90 * t87, t79, t80 * pkin(4) + t79 * qJ(5) + t81 * t90 + t82 * t88 + 0; t78, t88 * t87, t77, t78 * pkin(4) + t77 * qJ(5) + t81 * t88 - t82 * t90 + 0; t87 * t85, -t89, t87 * t84, -t89 * t91 + pkin(5) + 0 + (pkin(4) * t85 + qJ(5) * t84 + t83) * t87; 0, 0, 0, 1;];
	Tc_mdh = t1;
end