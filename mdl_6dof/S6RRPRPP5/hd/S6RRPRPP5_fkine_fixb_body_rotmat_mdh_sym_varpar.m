% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:07
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:07:38
	% EndTime: 2020-11-04 22:07:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:07:38
	% EndTime: 2020-11-04 22:07:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t49 = cos(qJ(1));
	t48 = sin(qJ(1));
	t1 = [t49, -t48, 0, 0; t48, t49, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:07:38
	% EndTime: 2020-11-04 22:07:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t53 = cos(qJ(1));
	t52 = cos(qJ(2));
	t51 = sin(qJ(1));
	t50 = sin(qJ(2));
	t1 = [t53 * t52, -t53 * t50, t51, t53 * pkin(1) + t51 * pkin(7) + 0; t51 * t52, -t51 * t50, -t53, t51 * pkin(1) - t53 * pkin(7) + 0; t50, t52, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:07:38
	% EndTime: 2020-11-04 22:07:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (16->14), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t58 = cos(qJ(1));
	t57 = cos(qJ(2));
	t56 = sin(qJ(1));
	t55 = sin(qJ(2));
	t54 = t57 * pkin(2) + t55 * qJ(3) + pkin(1);
	t1 = [t56, -t58 * t57, t58 * t55, t56 * pkin(7) + t54 * t58 + 0; -t58, -t56 * t57, t56 * t55, -t58 * pkin(7) + t54 * t56 + 0; 0, -t55, -t57, t55 * pkin(2) - t57 * qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:07:38
	% EndTime: 2020-11-04 22:07:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->14)
	t60 = sin(qJ(4));
	t62 = sin(qJ(1));
	t71 = t62 * t60;
	t63 = cos(qJ(4));
	t70 = t62 * t63;
	t65 = cos(qJ(1));
	t69 = t65 * t60;
	t68 = t65 * t63;
	t67 = pkin(2) + pkin(8);
	t66 = pkin(3) + pkin(7);
	t64 = cos(qJ(2));
	t61 = sin(qJ(2));
	t59 = t61 * qJ(3) + t67 * t64 + pkin(1);
	t1 = [t61 * t69 + t70, t61 * t68 - t71, t65 * t64, t59 * t65 + t66 * t62 + 0; t61 * t71 - t68, t61 * t70 + t69, t62 * t64, t59 * t62 - t66 * t65 + 0; -t64 * t60, -t64 * t63, t61, -t64 * qJ(3) + t67 * t61 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:07:38
	% EndTime: 2020-11-04 22:07:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->20), mult. (36->24), div. (0->0), fcn. (49->6), ass. (0->15)
	t75 = sin(qJ(4));
	t77 = sin(qJ(1));
	t85 = t77 * t75;
	t78 = cos(qJ(4));
	t84 = t77 * t78;
	t80 = cos(qJ(1));
	t83 = t80 * t75;
	t82 = t80 * t78;
	t81 = pkin(2) + pkin(8);
	t79 = cos(qJ(2));
	t76 = sin(qJ(2));
	t74 = -t75 * pkin(4) + t78 * qJ(5) - qJ(3);
	t73 = pkin(4) * t78 + qJ(5) * t75 + pkin(3) + pkin(7);
	t72 = -t74 * t76 + t81 * t79 + pkin(1);
	t1 = [t76 * t83 + t84, t80 * t79, -t76 * t82 + t85, t72 * t80 + t73 * t77 + 0; t76 * t85 - t82, t77 * t79, -t76 * t84 - t83, t72 * t77 - t73 * t80 + 0; -t79 * t75, t76, t79 * t78, t74 * t79 + t81 * t76 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:07:38
	% EndTime: 2020-11-04 22:07:38
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (42->24), mult. (36->24), div. (0->0), fcn. (49->6), ass. (0->16)
	t88 = sin(qJ(4));
	t90 = sin(qJ(1));
	t100 = t90 * t88;
	t91 = cos(qJ(4));
	t99 = t90 * t91;
	t93 = cos(qJ(1));
	t98 = t93 * t88;
	t97 = t93 * t91;
	t94 = pkin(4) + pkin(5);
	t96 = t91 * qJ(5) - t94 * t88 - qJ(3);
	t95 = qJ(5) * t88 + t94 * t91 + pkin(3) + pkin(7);
	t92 = cos(qJ(2));
	t89 = sin(qJ(2));
	t87 = -qJ(6) + pkin(2) + pkin(8);
	t86 = t87 * t92 - t96 * t89 + pkin(1);
	t1 = [t89 * t98 + t99, -t89 * t97 + t100, -t93 * t92, t86 * t93 + t95 * t90 + 0; t89 * t100 - t97, -t89 * t99 - t98, -t90 * t92, t86 * t90 - t95 * t93 + 0; -t92 * t88, t92 * t91, -t89, t87 * t89 + t96 * t92 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end