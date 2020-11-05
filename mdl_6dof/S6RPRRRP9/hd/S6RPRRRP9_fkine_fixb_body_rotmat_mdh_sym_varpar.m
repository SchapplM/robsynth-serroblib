% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRP9 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:53
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRRP9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:53:25
	% EndTime: 2020-11-04 21:53:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:53:25
	% EndTime: 2020-11-04 21:53:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t51 = cos(qJ(1));
	t50 = sin(qJ(1));
	t1 = [t51, -t50, 0, 0; t50, t51, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:53:25
	% EndTime: 2020-11-04 21:53:25
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t53 = cos(qJ(1));
	t52 = sin(qJ(1));
	t1 = [0, -t53, t52, pkin(1) * t53 + qJ(2) * t52 + 0; 0, -t52, -t53, pkin(1) * t52 - qJ(2) * t53 + 0; 1, 0, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:53:25
	% EndTime: 2020-11-04 21:53:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t58 = pkin(1) + pkin(7);
	t57 = cos(qJ(1));
	t56 = cos(qJ(3));
	t55 = sin(qJ(1));
	t54 = sin(qJ(3));
	t1 = [t55 * t54, t55 * t56, t57, t55 * qJ(2) + t58 * t57 + 0; -t57 * t54, -t57 * t56, t55, -t57 * qJ(2) + t58 * t55 + 0; t56, -t54, 0, pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:53:25
	% EndTime: 2020-11-04 21:53:25
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (20->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->13)
	t60 = sin(qJ(4));
	t62 = sin(qJ(1));
	t70 = t62 * t60;
	t63 = cos(qJ(4));
	t69 = t62 * t63;
	t65 = cos(qJ(1));
	t68 = t65 * t60;
	t67 = t65 * t63;
	t66 = pkin(1) + pkin(7);
	t64 = cos(qJ(3));
	t61 = sin(qJ(3));
	t59 = -t61 * pkin(3) + t64 * pkin(8) - qJ(2);
	t1 = [t61 * t69 + t68, -t61 * t70 + t67, -t62 * t64, -t59 * t62 + t66 * t65 + 0; -t61 * t67 + t70, t61 * t68 + t69, t65 * t64, t59 * t65 + t66 * t62 + 0; t64 * t63, -t64 * t60, t61, t64 * pkin(3) + t61 * pkin(8) + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:53:25
	% EndTime: 2020-11-04 21:53:25
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (38->21), mult. (31->22), div. (0->0), fcn. (44->8), ass. (0->16)
	t75 = qJ(4) + qJ(5);
	t73 = sin(t75);
	t77 = sin(qJ(1));
	t85 = t77 * t73;
	t74 = cos(t75);
	t84 = t77 * t74;
	t79 = cos(qJ(1));
	t83 = t79 * t73;
	t82 = t79 * t74;
	t72 = cos(qJ(4)) * pkin(4) + pkin(3);
	t76 = sin(qJ(3));
	t78 = cos(qJ(3));
	t80 = pkin(9) + pkin(8);
	t81 = t72 * t76 - t80 * t78 + qJ(2);
	t71 = sin(qJ(4)) * pkin(4) + pkin(1) + pkin(7);
	t1 = [t76 * t84 + t83, -t76 * t85 + t82, -t77 * t78, t71 * t79 + t81 * t77 + 0; -t76 * t82 + t85, t76 * t83 + t84, t79 * t78, t71 * t77 - t81 * t79 + 0; t78 * t74, -t78 * t73, t76, t78 * t72 + t80 * t76 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:53:25
	% EndTime: 2020-11-04 21:53:25
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (51->23), mult. (42->24), div. (0->0), fcn. (55->8), ass. (0->16)
	t91 = qJ(4) + qJ(5);
	t88 = sin(t91);
	t93 = sin(qJ(1));
	t102 = t93 * t88;
	t89 = cos(t91);
	t101 = t93 * t89;
	t95 = cos(qJ(1));
	t100 = t95 * t88;
	t99 = t95 * t89;
	t98 = pkin(5) * t88 + sin(qJ(4)) * pkin(4) + pkin(1) + pkin(7);
	t86 = pkin(5) * t89 + cos(qJ(4)) * pkin(4) + pkin(3);
	t90 = -qJ(6) - pkin(9) - pkin(8);
	t92 = sin(qJ(3));
	t94 = cos(qJ(3));
	t97 = t86 * t92 + t90 * t94 + qJ(2);
	t1 = [t92 * t101 + t100, -t92 * t102 + t99, -t93 * t94, t97 * t93 + t98 * t95 + 0; -t92 * t99 + t102, t92 * t100 + t101, t95 * t94, t98 * t93 - t97 * t95 + 0; t94 * t89, -t94 * t88, t92, t94 * t86 - t92 * t90 + pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end