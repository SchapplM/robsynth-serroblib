% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR6 (for one body)
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
% Datum: 2020-11-04 22:04
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:04:00
	% EndTime: 2020-11-04 22:04:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:04:00
	% EndTime: 2020-11-04 22:04:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t68 = cos(qJ(1));
	t67 = sin(qJ(1));
	t1 = [t68, -t67, 0, 0; t67, t68, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:04:00
	% EndTime: 2020-11-04 22:04:00
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t72 = cos(qJ(1));
	t71 = cos(qJ(2));
	t70 = sin(qJ(1));
	t69 = sin(qJ(2));
	t1 = [t72 * t71, -t72 * t69, t70, t72 * pkin(1) + t70 * pkin(7) + 0; t70 * t71, -t70 * t69, -t72, t70 * pkin(1) - t72 * pkin(7) + 0; t69, t71, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:04:00
	% EndTime: 2020-11-04 22:04:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (13->11), mult. (14->12), div. (0->0), fcn. (22->4), ass. (0->6)
	t77 = cos(qJ(1));
	t76 = cos(qJ(2));
	t75 = sin(qJ(1));
	t74 = sin(qJ(2));
	t73 = t76 * pkin(2) + t74 * qJ(3) + pkin(1);
	t1 = [t77 * t76, t75, t77 * t74, t75 * pkin(7) + t73 * t77 + 0; t75 * t76, -t77, t75 * t74, -t77 * pkin(7) + t73 * t75 + 0; t74, 0, -t76, t74 * pkin(2) - t76 * qJ(3) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:04:00
	% EndTime: 2020-11-04 22:04:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (23->16), mult. (26->16), div. (0->0), fcn. (40->6), ass. (0->12)
	t90 = cos(pkin(10));
	t89 = sin(pkin(10));
	t88 = pkin(2) + pkin(3);
	t87 = cos(qJ(1));
	t86 = cos(qJ(2));
	t85 = sin(qJ(1));
	t84 = sin(qJ(2));
	t83 = pkin(7) - qJ(4);
	t80 = t84 * qJ(3) + t88 * t86 + pkin(1);
	t79 = t84 * t90 - t86 * t89;
	t78 = -t84 * t89 - t86 * t90;
	t1 = [-t87 * t78, t87 * t79, -t85, t80 * t87 + t83 * t85 + 0; -t85 * t78, t85 * t79, t87, t80 * t85 - t83 * t87 + 0; t79, t78, 0, -t86 * qJ(3) + t88 * t84 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:04:00
	% EndTime: 2020-11-04 22:04:00
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (45->18), mult. (32->18), div. (0->0), fcn. (46->8), ass. (0->14)
	t99 = pkin(10) + qJ(5);
	t104 = sin(t99);
	t103 = cos(qJ(1));
	t102 = cos(qJ(2));
	t101 = sin(qJ(1));
	t100 = sin(qJ(2));
	t98 = qJ(4) - pkin(7) + pkin(8);
	t97 = cos(t99);
	t96 = sin(pkin(10)) * pkin(4) + qJ(3);
	t95 = cos(pkin(10)) * pkin(4) + pkin(2) + pkin(3);
	t93 = t100 * t97 - t102 * t104;
	t92 = t100 * t104 + t102 * t97;
	t91 = t96 * t100 + t95 * t102 + pkin(1);
	t1 = [t103 * t92, t103 * t93, -t101, -t98 * t101 + t91 * t103 + 0; t101 * t92, t101 * t93, t103, t91 * t101 + t98 * t103 + 0; t93, -t92, 0, t95 * t100 - t96 * t102 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:04:00
	% EndTime: 2020-11-04 22:04:00
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (77->34), mult. (75->37), div. (0->0), fcn. (97->15), ass. (0->26)
	t130 = pkin(2) + pkin(3);
	t116 = sin(qJ(6));
	t119 = sin(qJ(1));
	t129 = t119 * t116;
	t120 = cos(qJ(6));
	t128 = t119 * t120;
	t123 = cos(qJ(1));
	t127 = t123 * t116;
	t126 = t123 * t120;
	t125 = pkin(10) + qJ(5);
	t124 = sin(t125);
	t122 = cos(qJ(2));
	t121 = cos(qJ(5));
	t118 = sin(qJ(2));
	t117 = sin(qJ(5));
	t115 = cos(pkin(10));
	t114 = sin(pkin(10));
	t113 = qJ(4) - pkin(7) + pkin(8);
	t112 = -qJ(2) + t125;
	t111 = cos(t125);
	t109 = t115 * pkin(5) + t114 * pkin(9);
	t108 = -t114 * pkin(5) + t115 * pkin(9);
	t107 = t118 * t111 - t122 * t124;
	t106 = t122 * t111 + t118 * t124;
	t105 = (t115 * pkin(4) + t108 * t117 + t109 * t121 + t130) * t122 + (t114 * pkin(4) - t108 * t121 + t109 * t117 + qJ(3)) * t118 + pkin(1);
	t1 = [t106 * t126 - t129, -t106 * t127 - t128, -t123 * t107, t105 * t123 - t113 * t119 + 0; t106 * t128 + t127, -t106 * t129 + t126, -t119 * t107, t105 * t119 + t113 * t123 + 0; t107 * t120, -t107 * t116, t106, pkin(9) * cos(t112) - pkin(5) * sin(t112) - pkin(4) * sin(-qJ(2) + pkin(10)) + t130 * t118 - t122 * qJ(3) + 0 + pkin(6); 0, 0, 0, 1;];
	Tc_mdh = t1;
end