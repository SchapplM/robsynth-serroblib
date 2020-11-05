% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPP4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:07
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRPP4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPP4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:07:18
	% EndTime: 2020-11-04 22:07:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:07:18
	% EndTime: 2020-11-04 22:07:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t58 = cos(qJ(1));
	t57 = sin(qJ(1));
	t1 = [t58, -t57, 0, 0; t57, t58, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:07:18
	% EndTime: 2020-11-04 22:07:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t62 = cos(qJ(1));
	t61 = cos(qJ(2));
	t60 = sin(qJ(1));
	t59 = sin(qJ(2));
	t1 = [t62 * t61, -t62 * t59, t60, t62 * pkin(1) + t60 * pkin(7) + 0; t60 * t61, -t60 * t59, -t62, t60 * pkin(1) - t62 * pkin(7) + 0; t59, t61, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:07:18
	% EndTime: 2020-11-04 22:07:18
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (16->14), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t67 = cos(qJ(1));
	t66 = cos(qJ(2));
	t65 = sin(qJ(1));
	t64 = sin(qJ(2));
	t63 = t66 * pkin(2) + t64 * qJ(3) + pkin(1);
	t1 = [t65, -t67 * t66, t67 * t64, t65 * pkin(7) + t63 * t67 + 0; -t67, -t65 * t66, t65 * t64, -t67 * pkin(7) + t63 * t65 + 0; 0, -t64, -t66, t64 * pkin(2) - t66 * qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:07:18
	% EndTime: 2020-11-04 22:07:18
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->17), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->14)
	t69 = sin(qJ(4));
	t71 = sin(qJ(1));
	t80 = t71 * t69;
	t72 = cos(qJ(4));
	t79 = t71 * t72;
	t74 = cos(qJ(1));
	t78 = t74 * t69;
	t77 = t74 * t72;
	t76 = pkin(2) + pkin(8);
	t75 = pkin(3) + pkin(7);
	t73 = cos(qJ(2));
	t70 = sin(qJ(2));
	t68 = t70 * qJ(3) + t76 * t73 + pkin(1);
	t1 = [t70 * t78 + t79, t70 * t77 - t80, t74 * t73, t68 * t74 + t75 * t71 + 0; t70 * t80 - t77, t70 * t79 + t78, t71 * t73, t68 * t71 - t75 * t74 + 0; -t73 * t69, -t73 * t72, t70, -t73 * qJ(3) + t76 * t70 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:07:18
	% EndTime: 2020-11-04 22:07:18
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (40->20), mult. (31->22), div. (0->0), fcn. (44->8), ass. (0->16)
	t87 = qJ(4) + pkin(9);
	t84 = sin(t87);
	t89 = sin(qJ(1));
	t95 = t89 * t84;
	t85 = cos(t87);
	t94 = t89 * t85;
	t91 = cos(qJ(1));
	t93 = t91 * t84;
	t92 = t91 * t85;
	t90 = cos(qJ(2));
	t88 = sin(qJ(2));
	t86 = pkin(2) + pkin(8) + qJ(5);
	t83 = sin(qJ(4)) * pkin(4) + qJ(3);
	t82 = cos(qJ(4)) * pkin(4) + pkin(3) + pkin(7);
	t81 = t83 * t88 + t86 * t90 + pkin(1);
	t1 = [t88 * t93 + t94, t88 * t92 - t95, t91 * t90, t81 * t91 + t82 * t89 + 0; t88 * t95 - t92, t88 * t94 + t93, t89 * t90, t81 * t89 - t82 * t91 + 0; -t90 * t84, -t90 * t85, t88, -t83 * t90 + t86 * t88 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:07:18
	% EndTime: 2020-11-04 22:07:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (59->24), mult. (56->28), div. (0->0), fcn. (69->10), ass. (0->22)
	t107 = sin(qJ(1));
	t103 = qJ(4) + pkin(9);
	t99 = sin(t103);
	t117 = t107 * t99;
	t110 = cos(qJ(1));
	t116 = t110 * t99;
	t115 = sin(pkin(9));
	t100 = cos(t103);
	t114 = t107 * t100;
	t113 = t110 * t100;
	t105 = sin(qJ(4));
	t108 = cos(qJ(4));
	t104 = cos(pkin(9));
	t97 = pkin(5) * t104 + qJ(6) * t115 + pkin(4);
	t98 = -t115 * pkin(5) + qJ(6) * t104;
	t112 = t97 * t105 - t98 * t108 + qJ(3);
	t111 = t98 * t105 + t97 * t108 + pkin(3) + pkin(7);
	t109 = cos(qJ(2));
	t106 = sin(qJ(2));
	t102 = pkin(2) + pkin(8) + qJ(5);
	t96 = t102 * t109 + t112 * t106 + pkin(1);
	t1 = [t106 * t116 + t114, t110 * t109, -t106 * t113 + t117, t111 * t107 + t96 * t110 + 0; t106 * t117 - t113, t107 * t109, -t106 * t114 - t116, t96 * t107 - t111 * t110 + 0; -t109 * t99, t106, t109 * t100, t102 * t106 - t112 * t109 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end