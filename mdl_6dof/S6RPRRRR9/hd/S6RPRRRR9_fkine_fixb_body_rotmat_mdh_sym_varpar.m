% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRR9 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:56
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:56:53
	% EndTime: 2020-11-04 21:56:53
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:56:53
	% EndTime: 2020-11-04 21:56:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t50 = cos(qJ(1));
	t49 = sin(qJ(1));
	t1 = [t50, -t49, 0, 0; t49, t50, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:56:53
	% EndTime: 2020-11-04 21:56:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t52 = cos(qJ(1));
	t51 = sin(qJ(1));
	t1 = [0, -t52, t51, t52 * pkin(1) + t51 * qJ(2) + 0; 0, -t51, -t52, t51 * pkin(1) - t52 * qJ(2) + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:56:53
	% EndTime: 2020-11-04 21:56:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t57 = pkin(1) + pkin(7);
	t56 = cos(qJ(1));
	t55 = cos(qJ(3));
	t54 = sin(qJ(1));
	t53 = sin(qJ(3));
	t1 = [t54 * t53, t54 * t55, t56, t54 * qJ(2) + t57 * t56 + 0; -t56 * t53, -t56 * t55, t54, -t56 * qJ(2) + t57 * t54 + 0; t55, -t53, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:56:53
	% EndTime: 2020-11-04 21:56:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (20->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->13)
	t59 = sin(qJ(4));
	t61 = sin(qJ(1));
	t69 = t61 * t59;
	t62 = cos(qJ(4));
	t68 = t61 * t62;
	t64 = cos(qJ(1));
	t67 = t64 * t59;
	t66 = t64 * t62;
	t65 = pkin(1) + pkin(7);
	t63 = cos(qJ(3));
	t60 = sin(qJ(3));
	t58 = -t60 * pkin(3) + t63 * pkin(8) - qJ(2);
	t1 = [t60 * t68 + t67, -t60 * t69 + t66, -t61 * t63, -t58 * t61 + t65 * t64 + 0; -t60 * t66 + t69, t60 * t67 + t68, t64 * t63, t58 * t64 + t65 * t61 + 0; t63 * t62, -t63 * t59, t60, t63 * pkin(3) + t60 * pkin(8) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:56:53
	% EndTime: 2020-11-04 21:56:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (38->21), mult. (31->22), div. (0->0), fcn. (44->8), ass. (0->16)
	t74 = qJ(4) + qJ(5);
	t72 = sin(t74);
	t76 = sin(qJ(1));
	t84 = t76 * t72;
	t73 = cos(t74);
	t83 = t76 * t73;
	t78 = cos(qJ(1));
	t82 = t78 * t72;
	t81 = t78 * t73;
	t71 = cos(qJ(4)) * pkin(4) + pkin(3);
	t75 = sin(qJ(3));
	t77 = cos(qJ(3));
	t79 = pkin(9) + pkin(8);
	t80 = t71 * t75 - t79 * t77 + qJ(2);
	t70 = sin(qJ(4)) * pkin(4) + pkin(1) + pkin(7);
	t1 = [t75 * t83 + t82, -t75 * t84 + t81, -t76 * t77, t70 * t78 + t80 * t76 + 0; -t75 * t81 + t84, t75 * t82 + t83, t78 * t77, t70 * t76 - t80 * t78 + 0; t77 * t73, -t77 * t72, t75, t77 * t71 + t79 * t75 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:56:53
	% EndTime: 2020-11-04 21:56:53
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (61->24), mult. (42->24), div. (0->0), fcn. (55->10), ass. (0->17)
	t91 = qJ(4) + qJ(5);
	t89 = qJ(6) + t91;
	t87 = sin(t89);
	t93 = sin(qJ(1));
	t102 = t93 * t87;
	t88 = cos(t89);
	t101 = t93 * t88;
	t95 = cos(qJ(1));
	t100 = t95 * t87;
	t99 = t95 * t88;
	t98 = pkin(5) * sin(t91) + sin(qJ(4)) * pkin(4) + pkin(1) + pkin(7);
	t85 = pkin(5) * cos(t91) + cos(qJ(4)) * pkin(4) + pkin(3);
	t90 = -pkin(10) - pkin(9) - pkin(8);
	t92 = sin(qJ(3));
	t94 = cos(qJ(3));
	t97 = t85 * t92 + t90 * t94 + qJ(2);
	t1 = [t92 * t101 + t100, -t92 * t102 + t99, -t93 * t94, t97 * t93 + t98 * t95 + 0; -t92 * t99 + t102, t92 * t100 + t101, t95 * t94, t98 * t93 - t97 * t95 + 0; t94 * t88, -t94 * t87, t92, t94 * t85 - t92 * t90 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end