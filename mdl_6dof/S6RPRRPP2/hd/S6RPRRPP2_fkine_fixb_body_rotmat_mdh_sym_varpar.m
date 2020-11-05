% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPP2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:44
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRPP2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:44:15
	% EndTime: 2020-11-04 21:44:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:44:15
	% EndTime: 2020-11-04 21:44:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t61 = cos(qJ(1));
	t60 = sin(qJ(1));
	t1 = [t61, -t60, 0, 0; t60, t61, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:44:15
	% EndTime: 2020-11-04 21:44:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t64 = qJ(1) + pkin(9);
	t63 = cos(t64);
	t62 = sin(t64);
	t1 = [t63, -t62, 0, cos(qJ(1)) * pkin(1) + 0; t62, t63, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:44:15
	% EndTime: 2020-11-04 21:44:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t69 = cos(qJ(3));
	t68 = sin(qJ(3));
	t67 = qJ(1) + pkin(9);
	t66 = cos(t67);
	t65 = sin(t67);
	t1 = [t66 * t69, -t66 * t68, t65, t66 * pkin(2) + t65 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t65 * t69, -t65 * t68, -t66, t65 * pkin(2) - t66 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t68, t69, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:44:15
	% EndTime: 2020-11-04 21:44:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (38->19), mult. (32->24), div. (0->0), fcn. (45->8), ass. (0->11)
	t73 = sin(qJ(4));
	t76 = cos(qJ(3));
	t79 = t73 * t76;
	t75 = cos(qJ(4));
	t78 = t75 * t76;
	t74 = sin(qJ(3));
	t77 = pkin(3) * t76 + pkin(8) * t74 + pkin(2);
	t72 = qJ(1) + pkin(9);
	t71 = cos(t72);
	t70 = sin(t72);
	t1 = [t70 * t73 + t71 * t78, t70 * t75 - t71 * t79, t71 * t74, cos(qJ(1)) * pkin(1) + t70 * pkin(7) + 0 + t77 * t71; t70 * t78 - t71 * t73, -t70 * t79 - t71 * t75, t70 * t74, sin(qJ(1)) * pkin(1) - t71 * pkin(7) + 0 + t77 * t70; t74 * t75, -t74 * t73, -t76, t74 * pkin(3) - t76 * pkin(8) + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:44:15
	% EndTime: 2020-11-04 21:44:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (55->24), mult. (52->30), div. (0->0), fcn. (69->8), ass. (0->15)
	t87 = sin(qJ(4));
	t90 = cos(qJ(3));
	t93 = t87 * t90;
	t89 = cos(qJ(4));
	t92 = t89 * t90;
	t88 = sin(qJ(3));
	t91 = pkin(3) * t90 + pkin(8) * t88 + pkin(2);
	t86 = qJ(1) + pkin(9);
	t85 = cos(t86);
	t84 = sin(t86);
	t83 = t84 * t87 + t85 * t92;
	t82 = -t84 * t89 + t85 * t93;
	t81 = t84 * t92 - t85 * t87;
	t80 = t84 * t93 + t85 * t89;
	t1 = [t83, t85 * t88, t82, cos(qJ(1)) * pkin(1) + t83 * pkin(4) + t84 * pkin(7) + t82 * qJ(5) + 0 + t91 * t85; t81, t84 * t88, t80, sin(qJ(1)) * pkin(1) + t81 * pkin(4) - t85 * pkin(7) + t80 * qJ(5) + 0 + t91 * t84; t88 * t89, -t90, t88 * t87, -t90 * pkin(8) + pkin(6) + qJ(2) + 0 + (pkin(4) * t89 + qJ(5) * t87 + pkin(3)) * t88; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:44:15
	% EndTime: 2020-11-04 21:44:15
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (70->26), mult. (72->30), div. (0->0), fcn. (85->10), ass. (0->21)
	t103 = sin(qJ(4));
	t107 = cos(qJ(3));
	t113 = t103 * t107;
	t106 = cos(qJ(4));
	t112 = t106 * t107;
	t102 = qJ(6) - pkin(8);
	t104 = sin(qJ(3));
	t109 = pkin(4) + pkin(5);
	t96 = t103 * qJ(5) + t109 * t106 + pkin(3);
	t111 = -t102 * t104 + t107 * t96 + pkin(2);
	t110 = qJ(5) * t106 - t103 * t109 - pkin(7);
	t108 = cos(qJ(1));
	t105 = sin(qJ(1));
	t101 = cos(pkin(9));
	t100 = sin(pkin(9));
	t99 = qJ(1) + pkin(9);
	t98 = cos(t99);
	t97 = sin(t99);
	t95 = t111 * t100 + t110 * t101;
	t94 = -t110 * t100 + t111 * t101 + pkin(1);
	t1 = [t97 * t103 + t98 * t112, -t97 * t106 + t98 * t113, -t98 * t104, -t95 * t105 + t94 * t108 + 0; -t98 * t103 + t97 * t112, t98 * t106 + t97 * t113, -t97 * t104, t94 * t105 + t95 * t108 + 0; t104 * t106, t104 * t103, t107, t102 * t107 + t96 * t104 + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end