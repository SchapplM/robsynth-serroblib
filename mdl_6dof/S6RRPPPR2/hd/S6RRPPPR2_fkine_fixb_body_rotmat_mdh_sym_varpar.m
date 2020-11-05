% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPPR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:59
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPPR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:03
	% EndTime: 2020-11-04 21:59:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:03
	% EndTime: 2020-11-04 21:59:03
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t1 = [t57, -t56, 0, 0; t56, t57, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:03
	% EndTime: 2020-11-04 21:59:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t61 = cos(qJ(1));
	t60 = cos(qJ(2));
	t59 = sin(qJ(1));
	t58 = sin(qJ(2));
	t1 = [t61 * t60, -t61 * t58, t59, t61 * pkin(1) + t59 * pkin(7) + 0; t59 * t60, -t59 * t58, -t61, t59 * pkin(1) - t61 * pkin(7) + 0; t58, t60, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:03
	% EndTime: 2020-11-04 21:59:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t68 = cos(qJ(1));
	t67 = sin(qJ(1));
	t66 = -qJ(3) - pkin(7);
	t65 = qJ(2) + pkin(9);
	t64 = cos(t65);
	t63 = sin(t65);
	t62 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t68 * t64, -t68 * t63, t67, t68 * t62 - t66 * t67 + 0; t67 * t64, -t67 * t63, -t68, t67 * t62 + t68 * t66 + 0; t63, t64, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:03
	% EndTime: 2020-11-04 21:59:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->20), mult. (23->17), div. (0->0), fcn. (31->8), ass. (0->11)
	t78 = cos(qJ(1));
	t77 = sin(qJ(1));
	t76 = sin(qJ(2));
	t75 = -qJ(3) - pkin(7);
	t74 = cos(pkin(9));
	t73 = sin(pkin(9));
	t72 = qJ(2) + pkin(9);
	t71 = cos(t72);
	t70 = sin(t72);
	t69 = (pkin(3) * t74 + qJ(4) * t73 + pkin(2)) * cos(qJ(2)) + (-t73 * pkin(3) + qJ(4) * t74) * t76 + pkin(1);
	t1 = [t77, -t78 * t71, t78 * t70, t69 * t78 - t75 * t77 + 0; -t78, -t77 * t71, t77 * t70, t69 * t77 + t78 * t75 + 0; 0, -t70, -t71, t76 * pkin(2) + t70 * pkin(3) - t71 * qJ(4) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:03
	% EndTime: 2020-11-04 21:59:03
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (44->22), mult. (35->25), div. (0->0), fcn. (48->10), ass. (0->18)
	t84 = sin(pkin(10));
	t90 = sin(qJ(1));
	t95 = t90 * t84;
	t86 = cos(pkin(10));
	t94 = t90 * t86;
	t91 = cos(qJ(1));
	t93 = t91 * t84;
	t92 = t91 * t86;
	t89 = sin(qJ(2));
	t88 = pkin(3) + qJ(5);
	t87 = cos(pkin(9));
	t85 = sin(pkin(9));
	t83 = qJ(2) + pkin(9);
	t82 = qJ(3) + pkin(4) + pkin(7);
	t81 = cos(t83);
	t80 = sin(t83);
	t79 = (qJ(4) * t85 + t88 * t87 + pkin(2)) * cos(qJ(2)) + (qJ(4) * t87 - t85 * t88) * t89 + pkin(1);
	t1 = [t80 * t93 + t94, t80 * t92 - t95, t91 * t81, t79 * t91 + t82 * t90 + 0; t80 * t95 - t92, t80 * t94 + t93, t90 * t81, t79 * t90 - t82 * t91 + 0; -t81 * t84, -t81 * t86, t80, t89 * pkin(2) - t81 * qJ(4) + t88 * t80 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:59:03
	% EndTime: 2020-11-04 21:59:03
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (71->30), mult. (49->30), div. (0->0), fcn. (56->14), ass. (0->20)
	t109 = sin(qJ(1));
	t104 = pkin(10) + qJ(6);
	t99 = sin(t104);
	t114 = t109 * t99;
	t110 = cos(qJ(1));
	t113 = t110 * t99;
	t101 = cos(t104);
	t112 = t109 * t101;
	t111 = t110 * t101;
	t105 = qJ(2) + pkin(9);
	t108 = sin(qJ(2));
	t107 = cos(pkin(9));
	t106 = sin(pkin(9));
	t103 = qJ(5) + pkin(3) + pkin(8);
	t102 = cos(t105);
	t100 = sin(t105);
	t98 = sin(pkin(10)) * pkin(5) + qJ(4);
	t97 = cos(pkin(10)) * pkin(5) + qJ(3) + pkin(4) + pkin(7);
	t96 = (t103 * t107 + t98 * t106 + pkin(2)) * cos(qJ(2)) + (-t106 * t103 + t98 * t107) * t108 + pkin(1);
	t1 = [t100 * t113 + t112, t100 * t111 - t114, t110 * t102, t97 * t109 + t96 * t110 + 0; t100 * t114 - t111, t100 * t112 + t113, t109 * t102, t96 * t109 - t97 * t110 + 0; -t102 * t99, -t102 * t101, t100, t103 * t100 - t102 * qJ(4) + t108 * pkin(2) + 0 + pkin(6) + (sin(-pkin(10) + t105) / 0.2e1 - sin(pkin(10) + t105) / 0.2e1) * pkin(5); 0, 0, 0, 1;];
	Tc_mdh = t1;
end