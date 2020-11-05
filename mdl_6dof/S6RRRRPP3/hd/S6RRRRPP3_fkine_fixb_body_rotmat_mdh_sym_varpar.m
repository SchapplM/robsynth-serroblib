% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPP3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:35
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:35:39
	% EndTime: 2020-11-04 22:35:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:35:39
	% EndTime: 2020-11-04 22:35:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t64 = cos(qJ(1));
	t63 = sin(qJ(1));
	t1 = [t64, -t63, 0, 0; t63, t64, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:35:39
	% EndTime: 2020-11-04 22:35:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t68 = cos(qJ(1));
	t67 = cos(qJ(2));
	t66 = sin(qJ(1));
	t65 = sin(qJ(2));
	t1 = [t68 * t67, -t68 * t65, t66, t68 * pkin(1) + t66 * pkin(7) + 0; t66 * t67, -t66 * t65, -t68, t66 * pkin(1) - t68 * pkin(7) + 0; t65, t67, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:35:39
	% EndTime: 2020-11-04 22:35:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t75 = pkin(8) + pkin(7);
	t74 = cos(qJ(1));
	t73 = sin(qJ(1));
	t72 = qJ(2) + qJ(3);
	t71 = cos(t72);
	t70 = sin(t72);
	t69 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t74 * t71, -t74 * t70, t73, t74 * t69 + t75 * t73 + 0; t73 * t71, -t73 * t70, -t74, t73 * t69 - t74 * t75 + 0; t70, t71, 0, sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:35:39
	% EndTime: 2020-11-04 22:35:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t80 = sin(qJ(4));
	t81 = sin(qJ(1));
	t89 = t81 * t80;
	t82 = cos(qJ(4));
	t88 = t81 * t82;
	t83 = cos(qJ(1));
	t87 = t83 * t80;
	t86 = t83 * t82;
	t79 = qJ(2) + qJ(3);
	t77 = sin(t79);
	t78 = cos(t79);
	t85 = pkin(3) * t78 + pkin(9) * t77 + cos(qJ(2)) * pkin(2) + pkin(1);
	t84 = pkin(8) + pkin(7);
	t1 = [t78 * t86 + t89, -t78 * t87 + t88, t83 * t77, t84 * t81 + t85 * t83 + 0; t78 * t88 - t87, -t78 * t89 - t86, t81 * t77, t85 * t81 - t83 * t84 + 0; t77 * t82, -t77 * t80, -t78, t77 * pkin(3) - t78 * pkin(9) + sin(qJ(2)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:35:39
	% EndTime: 2020-11-04 22:35:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (53->27), mult. (53->28), div. (0->0), fcn. (70->8), ass. (0->18)
	t98 = sin(qJ(4));
	t99 = sin(qJ(1));
	t107 = t99 * t98;
	t101 = cos(qJ(1));
	t106 = t101 * t98;
	t100 = cos(qJ(4));
	t105 = t99 * t100;
	t104 = t101 * t100;
	t97 = qJ(2) + qJ(3);
	t95 = sin(t97);
	t96 = cos(t97);
	t103 = pkin(3) * t96 + pkin(9) * t95 + cos(qJ(2)) * pkin(2) + pkin(1);
	t102 = pkin(8) + pkin(7);
	t93 = t96 * t104 + t107;
	t92 = t96 * t106 - t105;
	t91 = t96 * t105 - t106;
	t90 = t96 * t107 + t104;
	t1 = [t101 * t95, -t93, t92, t93 * pkin(4) + t92 * qJ(5) + t103 * t101 + t102 * t99 + 0; t99 * t95, -t91, t90, t91 * pkin(4) + t90 * qJ(5) - t101 * t102 + t103 * t99 + 0; -t96, -t95 * t100, t95 * t98, sin(qJ(2)) * pkin(2) - t96 * pkin(9) + pkin(6) + 0 + (pkin(4) * t100 + qJ(5) * t98 + pkin(3)) * t95; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:35:39
	% EndTime: 2020-11-04 22:35:39
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (73->32), mult. (61->35), div. (0->0), fcn. (67->14), ass. (0->22)
	t116 = sin(qJ(4));
	t119 = sin(qJ(1));
	t128 = t119 * t116;
	t120 = cos(qJ(4));
	t127 = t119 * t120;
	t122 = cos(qJ(1));
	t126 = t122 * t116;
	t125 = t122 * t120;
	t114 = qJ(2) + qJ(3);
	t115 = qJ(6) + pkin(4);
	t124 = qJ(5) * t120 - t115 * t116 - pkin(7) - pkin(8);
	t123 = pkin(5) + pkin(9);
	t121 = cos(qJ(3));
	t118 = sin(qJ(2));
	t117 = sin(qJ(3));
	t113 = -qJ(4) + t114;
	t112 = qJ(4) + t114;
	t111 = cos(t114);
	t110 = sin(t114);
	t109 = qJ(5) * t116 + t115 * t120 + pkin(3);
	t108 = (t109 * t121 + t123 * t117 + pkin(2)) * cos(qJ(2)) + pkin(1) + (-t109 * t117 + t123 * t121) * t118;
	t1 = [t122 * t110, t111 * t126 - t127, t111 * t125 + t128, t108 * t122 - t124 * t119 + 0; t119 * t110, t111 * t128 + t125, t111 * t127 - t126, t108 * t119 + t124 * t122 + 0; -t111, t110 * t116, t110 * t120, -t123 * t111 + t118 * pkin(2) + t110 * pkin(3) + 0 + pkin(6) + (sin(t113) / 0.2e1 + sin(t112) / 0.2e1) * t115 + (cos(t113) / 0.2e1 - cos(t112) / 0.2e1) * qJ(5); 0, 0, 0, 1;];
	Tc_mdh = t1;
end