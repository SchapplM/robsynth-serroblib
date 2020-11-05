% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRP2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:00
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPRP2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:00:47
	% EndTime: 2020-11-04 22:00:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:00:47
	% EndTime: 2020-11-04 22:00:47
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t1 = [t56, -t55, 0, 0; t55, t56, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:00:47
	% EndTime: 2020-11-04 22:00:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t60 = cos(qJ(1));
	t59 = cos(qJ(2));
	t58 = sin(qJ(1));
	t57 = sin(qJ(2));
	t1 = [t60 * t59, -t60 * t57, t58, t60 * pkin(1) + t58 * pkin(7) + 0; t58 * t59, -t58 * t57, -t60, t58 * pkin(1) - t60 * pkin(7) + 0; t57, t59, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:00:47
	% EndTime: 2020-11-04 22:00:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t67 = cos(qJ(1));
	t66 = sin(qJ(1));
	t65 = -qJ(3) - pkin(7);
	t64 = qJ(2) + pkin(9);
	t63 = cos(t64);
	t62 = sin(t64);
	t61 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t67 * t63, -t67 * t62, t66, t67 * t61 - t65 * t66 + 0; t66 * t63, -t66 * t62, -t67, t66 * t61 + t67 * t65 + 0; t62, t63, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:00:47
	% EndTime: 2020-11-04 22:00:47
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (33->20), mult. (23->17), div. (0->0), fcn. (31->8), ass. (0->11)
	t77 = cos(qJ(1));
	t76 = sin(qJ(1));
	t75 = sin(qJ(2));
	t74 = -qJ(3) - pkin(7);
	t73 = cos(pkin(9));
	t72 = sin(pkin(9));
	t71 = qJ(2) + pkin(9);
	t70 = cos(t71);
	t69 = sin(t71);
	t68 = (pkin(3) * t73 + qJ(4) * t72 + pkin(2)) * cos(qJ(2)) + (-t72 * pkin(3) + qJ(4) * t73) * t75 + pkin(1);
	t1 = [t76, -t77 * t70, t77 * t69, t68 * t77 - t74 * t76 + 0; -t77, -t76 * t70, t76 * t69, t68 * t76 + t77 * t74 + 0; 0, -t69, -t70, t75 * pkin(2) + t69 * pkin(3) - t70 * qJ(4) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:00:47
	% EndTime: 2020-11-04 22:00:47
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (44->22), mult. (35->25), div. (0->0), fcn. (48->10), ass. (0->18)
	t85 = sin(qJ(5));
	t87 = sin(qJ(1));
	t94 = t87 * t85;
	t88 = cos(qJ(5));
	t93 = t87 * t88;
	t89 = cos(qJ(1));
	t92 = t89 * t85;
	t91 = t89 * t88;
	t90 = pkin(3) + pkin(8);
	t86 = sin(qJ(2));
	t84 = cos(pkin(9));
	t83 = sin(pkin(9));
	t82 = qJ(2) + pkin(9);
	t81 = qJ(3) + pkin(4) + pkin(7);
	t80 = cos(t82);
	t79 = sin(t82);
	t78 = (qJ(4) * t83 + t90 * t84 + pkin(2)) * cos(qJ(2)) + (qJ(4) * t84 - t83 * t90) * t86 + pkin(1);
	t1 = [t79 * t92 + t93, t79 * t91 - t94, t89 * t80, t78 * t89 + t81 * t87 + 0; t79 * t94 - t91, t79 * t93 + t92, t87 * t80, t78 * t87 - t81 * t89 + 0; -t80 * t85, -t80 * t88, t79, t86 * pkin(2) - t80 * qJ(4) + t90 * t79 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:00:47
	% EndTime: 2020-11-04 22:00:47
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (61->29), mult. (53->30), div. (0->0), fcn. (60->12), ass. (0->19)
	t103 = sin(qJ(5));
	t105 = sin(qJ(1));
	t112 = t105 * t103;
	t106 = cos(qJ(5));
	t111 = t105 * t106;
	t107 = cos(qJ(1));
	t110 = t107 * t103;
	t109 = t107 * t106;
	t100 = qJ(2) + pkin(9);
	t108 = pkin(5) * t103 + qJ(4);
	t104 = sin(qJ(2));
	t102 = cos(pkin(9));
	t101 = sin(pkin(9));
	t99 = qJ(6) + pkin(3) + pkin(8);
	t98 = cos(t100);
	t97 = sin(t100);
	t96 = t106 * pkin(5) + pkin(4) + pkin(7) + qJ(3);
	t95 = (t108 * t101 + t99 * t102 + pkin(2)) * cos(qJ(2)) + (-t101 * t99 + t108 * t102) * t104 + pkin(1);
	t1 = [t97 * t110 + t111, t97 * t109 - t112, t107 * t98, t96 * t105 + t95 * t107 + 0; t97 * t112 - t109, t97 * t111 + t110, t105 * t98, t95 * t105 - t96 * t107 + 0; -t98 * t103, -t98 * t106, t97, t99 * t97 + t104 * pkin(2) - t98 * qJ(4) + 0 + pkin(6) + (sin(-qJ(5) + t100) / 0.2e1 - sin(qJ(5) + t100) / 0.2e1) * pkin(5); 0, 0, 0, 1;];
	Tc_mdh = t1;
end