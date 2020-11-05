% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:05
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPRR10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:32
	% EndTime: 2020-11-04 22:05:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:32
	% EndTime: 2020-11-04 22:05:32
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t51 = cos(qJ(1));
	t50 = sin(qJ(1));
	t1 = [t51, -t50, 0, 0; t50, t51, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:32
	% EndTime: 2020-11-04 22:05:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t55 = cos(qJ(1));
	t54 = cos(qJ(2));
	t53 = sin(qJ(1));
	t52 = sin(qJ(2));
	t1 = [t55 * t54, -t55 * t52, t53, t55 * pkin(1) + t53 * pkin(7) + 0; t53 * t54, -t53 * t52, -t55, t53 * pkin(1) - t55 * pkin(7) + 0; t52, t54, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:32
	% EndTime: 2020-11-04 22:05:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (16->14), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t60 = cos(qJ(1));
	t59 = cos(qJ(2));
	t58 = sin(qJ(1));
	t57 = sin(qJ(2));
	t56 = t59 * pkin(2) + t57 * qJ(3) + pkin(1);
	t1 = [t58, -t60 * t59, t60 * t57, t58 * pkin(7) + t56 * t60 + 0; -t60, -t58 * t59, t58 * t57, -t60 * pkin(7) + t56 * t58 + 0; 0, -t57, -t59, t57 * pkin(2) - t59 * qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:32
	% EndTime: 2020-11-04 22:05:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->14)
	t62 = sin(pkin(10));
	t66 = sin(qJ(1));
	t73 = t66 * t62;
	t63 = cos(pkin(10));
	t72 = t66 * t63;
	t68 = cos(qJ(1));
	t71 = t68 * t62;
	t70 = t68 * t63;
	t69 = pkin(3) + pkin(7);
	t67 = cos(qJ(2));
	t65 = sin(qJ(2));
	t64 = pkin(2) + qJ(4);
	t61 = t65 * qJ(3) + t64 * t67 + pkin(1);
	t1 = [t65 * t71 + t72, t65 * t70 - t73, t68 * t67, t61 * t68 + t69 * t66 + 0; t65 * t73 - t70, t65 * t72 + t71, t66 * t67, t61 * t66 - t69 * t68 + 0; -t67 * t62, -t67 * t63, t65, -t67 * qJ(3) + t64 * t65 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:32
	% EndTime: 2020-11-04 22:05:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (40->20), mult. (31->22), div. (0->0), fcn. (44->8), ass. (0->16)
	t80 = pkin(10) + qJ(5);
	t77 = sin(t80);
	t82 = sin(qJ(1));
	t88 = t82 * t77;
	t78 = cos(t80);
	t87 = t82 * t78;
	t84 = cos(qJ(1));
	t86 = t84 * t77;
	t85 = t84 * t78;
	t83 = cos(qJ(2));
	t81 = sin(qJ(2));
	t79 = pkin(2) + pkin(8) + qJ(4);
	t76 = sin(pkin(10)) * pkin(4) + qJ(3);
	t75 = cos(pkin(10)) * pkin(4) + pkin(3) + pkin(7);
	t74 = t76 * t81 + t79 * t83 + pkin(1);
	t1 = [t81 * t86 + t87, t81 * t85 - t88, t84 * t83, t74 * t84 + t75 * t82 + 0; t81 * t88 - t85, t81 * t87 + t86, t82 * t83, t74 * t82 - t75 * t84 + 0; -t83 * t77, -t83 * t78, t81, -t76 * t83 + t79 * t81 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:05:32
	% EndTime: 2020-11-04 22:05:32
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (63->26), mult. (49->30), div. (0->0), fcn. (62->11), ass. (0->20)
	t104 = pkin(10) + qJ(5);
	t93 = qJ(6) + t104;
	t91 = sin(t93);
	t99 = sin(qJ(1));
	t107 = t99 * t91;
	t92 = cos(t93);
	t106 = t99 * t92;
	t102 = cos(qJ(1));
	t98 = sin(qJ(2));
	t105 = t102 * t98;
	t95 = sin(pkin(10));
	t103 = t95 * pkin(4) + qJ(3);
	t101 = cos(qJ(2));
	t100 = cos(qJ(5));
	t97 = sin(qJ(5));
	t96 = cos(pkin(10));
	t94 = qJ(4) + pkin(2) + pkin(8) + pkin(9);
	t90 = t96 * pkin(4) + pkin(3) + pkin(7) + (t100 * t96 - t95 * t97) * pkin(5);
	t89 = t94 * t101 + ((t100 * t95 + t96 * t97) * pkin(5) + t103) * t98 + pkin(1);
	t1 = [t91 * t105 + t106, t92 * t105 - t107, t102 * t101, t89 * t102 + t90 * t99 + 0; -t102 * t92 + t98 * t107, t102 * t91 + t98 * t106, t99 * t101, -t90 * t102 + t89 * t99 + 0; -t101 * t91, -t101 * t92, t98, t94 * t98 + pkin(6) + 0 + (-sin(t104) * pkin(5) - t103) * t101; 0, 0, 0, 1;];
	Tc_mdh = t1;
end