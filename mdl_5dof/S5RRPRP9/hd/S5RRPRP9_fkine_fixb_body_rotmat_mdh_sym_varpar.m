% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRP9 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:34
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPRP9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:34:40
	% EndTime: 2020-11-04 20:34:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:34:40
	% EndTime: 2020-11-04 20:34:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t48 = cos(qJ(1));
	t47 = sin(qJ(1));
	t1 = [t48, -t47, 0, 0; t47, t48, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:34:40
	% EndTime: 2020-11-04 20:34:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t52 = cos(qJ(1));
	t51 = cos(qJ(2));
	t50 = sin(qJ(1));
	t49 = sin(qJ(2));
	t1 = [t52 * t51, -t52 * t49, t50, t52 * pkin(1) + t50 * pkin(6) + 0; t50 * t51, -t50 * t49, -t52, t50 * pkin(1) - t52 * pkin(6) + 0; t49, t51, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:34:40
	% EndTime: 2020-11-04 20:34:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (17->15), mult. (26->21), div. (0->0), fcn. (39->6), ass. (0->11)
	t57 = sin(qJ(1));
	t58 = cos(qJ(2));
	t62 = t57 * t58;
	t54 = sin(pkin(8));
	t59 = cos(qJ(1));
	t61 = t59 * t54;
	t55 = cos(pkin(8));
	t60 = t59 * t55;
	t56 = sin(qJ(2));
	t53 = t58 * pkin(2) + t56 * qJ(3) + pkin(1);
	t1 = [t57 * t54 + t58 * t60, t57 * t55 - t58 * t61, t59 * t56, t57 * pkin(6) + t53 * t59 + 0; t55 * t62 - t61, -t54 * t62 - t60, t57 * t56, -t59 * pkin(6) + t53 * t57 + 0; t56 * t55, -t56 * t54, -t58, t56 * pkin(2) - t58 * qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:34:40
	% EndTime: 2020-11-04 20:34:40
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (35->19), mult. (31->23), div. (0->0), fcn. (44->8), ass. (0->15)
	t71 = sin(qJ(1));
	t72 = cos(qJ(2));
	t76 = t71 * t72;
	t68 = pkin(8) + qJ(4);
	t66 = sin(t68);
	t73 = cos(qJ(1));
	t75 = t73 * t66;
	t67 = cos(t68);
	t74 = t73 * t67;
	t70 = sin(qJ(2));
	t69 = qJ(3) + pkin(7);
	t65 = cos(pkin(8)) * pkin(3) + pkin(2);
	t64 = sin(pkin(8)) * pkin(3) + pkin(6);
	t63 = t65 * t72 + t69 * t70 + pkin(1);
	t1 = [t71 * t66 + t72 * t74, t71 * t67 - t72 * t75, t73 * t70, t63 * t73 + t64 * t71 + 0; t67 * t76 - t75, -t66 * t76 - t74, t71 * t70, t63 * t71 - t64 * t73 + 0; t70 * t67, -t70 * t66, -t72, t70 * t65 - t72 * t69 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:34:40
	% EndTime: 2020-11-04 20:34:40
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (54->24), mult. (51->29), div. (0->0), fcn. (68->8), ass. (0->19)
	t89 = sin(qJ(1));
	t90 = cos(qJ(2));
	t94 = t89 * t90;
	t86 = pkin(8) + qJ(4);
	t84 = sin(t86);
	t91 = cos(qJ(1));
	t93 = t91 * t84;
	t85 = cos(t86);
	t92 = t91 * t85;
	t88 = sin(qJ(2));
	t87 = qJ(3) + pkin(7);
	t83 = cos(pkin(8)) * pkin(3) + pkin(2);
	t82 = sin(pkin(8)) * pkin(3) + pkin(6);
	t81 = t83 * t90 + t87 * t88 + pkin(1);
	t80 = t89 * t84 + t90 * t92;
	t79 = -t89 * t85 + t90 * t93;
	t78 = t85 * t94 - t93;
	t77 = t84 * t94 + t92;
	t1 = [t80, t91 * t88, t79, t80 * pkin(4) + t79 * qJ(5) + t81 * t91 + t82 * t89 + 0; t78, t89 * t88, t77, t78 * pkin(4) + t77 * qJ(5) + t81 * t89 - t82 * t91 + 0; t88 * t85, -t90, t88 * t84, -t90 * t87 + pkin(5) + 0 + (pkin(4) * t85 + qJ(5) * t84 + t83) * t88; 0, 0, 0, 1;];
	Tc_mdh = t1;
end