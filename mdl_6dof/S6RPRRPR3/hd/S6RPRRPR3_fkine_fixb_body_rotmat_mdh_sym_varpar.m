% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:47
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:47:17
	% EndTime: 2020-11-04 21:47:17
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:47:17
	% EndTime: 2020-11-04 21:47:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t67 = cos(qJ(1));
	t66 = sin(qJ(1));
	t1 = [t67, -t66, 0, 0; t66, t67, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:47:17
	% EndTime: 2020-11-04 21:47:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t70 = qJ(1) + pkin(10);
	t69 = cos(t70);
	t68 = sin(t70);
	t1 = [t69, -t68, 0, cos(qJ(1)) * pkin(1) + 0; t68, t69, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:47:17
	% EndTime: 2020-11-04 21:47:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t75 = cos(qJ(3));
	t74 = sin(qJ(3));
	t73 = qJ(1) + pkin(10);
	t72 = cos(t73);
	t71 = sin(t73);
	t1 = [t72 * t75, -t72 * t74, t71, t72 * pkin(2) + t71 * pkin(7) + cos(qJ(1)) * pkin(1) + 0; t71 * t75, -t71 * t74, -t72, t71 * pkin(2) - t72 * pkin(7) + sin(qJ(1)) * pkin(1) + 0; t74, t75, 0, qJ(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:47:17
	% EndTime: 2020-11-04 21:47:17
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (38->19), mult. (32->24), div. (0->0), fcn. (45->8), ass. (0->11)
	t79 = sin(qJ(4));
	t82 = cos(qJ(3));
	t85 = t79 * t82;
	t81 = cos(qJ(4));
	t84 = t81 * t82;
	t80 = sin(qJ(3));
	t83 = pkin(3) * t82 + pkin(8) * t80 + pkin(2);
	t78 = qJ(1) + pkin(10);
	t77 = cos(t78);
	t76 = sin(t78);
	t1 = [t76 * t79 + t77 * t84, t76 * t81 - t77 * t85, t77 * t80, cos(qJ(1)) * pkin(1) + t76 * pkin(7) + 0 + t83 * t77; t76 * t84 - t77 * t79, -t76 * t85 - t77 * t81, t76 * t80, sin(qJ(1)) * pkin(1) - t77 * pkin(7) + 0 + t83 * t76; t80 * t81, -t80 * t79, -t82, t80 * pkin(3) - t82 * pkin(8) + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:47:17
	% EndTime: 2020-11-04 21:47:17
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (55->24), mult. (52->30), div. (0->0), fcn. (69->8), ass. (0->15)
	t93 = sin(qJ(4));
	t96 = cos(qJ(3));
	t99 = t93 * t96;
	t95 = cos(qJ(4));
	t98 = t95 * t96;
	t94 = sin(qJ(3));
	t97 = pkin(3) * t96 + pkin(8) * t94 + pkin(2);
	t92 = qJ(1) + pkin(10);
	t91 = cos(t92);
	t90 = sin(t92);
	t89 = t90 * t93 + t91 * t98;
	t88 = -t90 * t95 + t91 * t99;
	t87 = t90 * t98 - t91 * t93;
	t86 = t90 * t99 + t91 * t95;
	t1 = [t89, t91 * t94, t88, cos(qJ(1)) * pkin(1) + t89 * pkin(4) + t90 * pkin(7) + t88 * qJ(5) + 0 + t97 * t91; t87, t90 * t94, t86, sin(qJ(1)) * pkin(1) + t87 * pkin(4) - t91 * pkin(7) + t86 * qJ(5) + 0 + t97 * t90; t94 * t95, -t96, t94 * t93, -t96 * pkin(8) + pkin(6) + qJ(2) + 0 + (pkin(4) * t95 + qJ(5) * t93 + pkin(3)) * t94; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:47:17
	% EndTime: 2020-11-04 21:47:17
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (80->28), mult. (92->34), div. (0->0), fcn. (115->12), ass. (0->25)
	t107 = qJ(1) + pkin(10);
	t105 = sin(t107);
	t116 = cos(qJ(3));
	t123 = t105 * t116;
	t106 = cos(t107);
	t122 = t106 * t116;
	t111 = sin(qJ(4));
	t115 = cos(qJ(4));
	t119 = pkin(4) + pkin(5);
	t104 = t111 * qJ(5) + t119 * t115 + pkin(3);
	t112 = sin(qJ(3));
	t118 = pkin(8) - pkin(9);
	t121 = t104 * t116 + t112 * t118 + pkin(2);
	t120 = qJ(5) * t115 - t111 * t119 - pkin(7);
	t117 = cos(qJ(1));
	t114 = cos(qJ(6));
	t113 = sin(qJ(1));
	t110 = sin(qJ(6));
	t109 = cos(pkin(10));
	t108 = sin(pkin(10));
	t103 = t110 * t111 + t114 * t115;
	t102 = -t110 * t115 + t114 * t111;
	t101 = t121 * t108 + t120 * t109;
	t100 = -t120 * t108 + t121 * t109 + pkin(1);
	t1 = [t105 * t102 + t103 * t122, t102 * t122 - t105 * t103, -t106 * t112, t100 * t117 - t101 * t113 + 0; -t102 * t106 + t103 * t123, t102 * t123 + t103 * t106, -t105 * t112, t100 * t113 + t101 * t117 + 0; t112 * t103, t112 * t102, t116, t104 * t112 - t118 * t116 + pkin(6) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end