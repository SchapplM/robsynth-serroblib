% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPPR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:34
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:06
	% EndTime: 2020-11-04 21:34:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:06
	% EndTime: 2020-11-04 21:34:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t70 = cos(qJ(1));
	t69 = sin(qJ(1));
	t1 = [t70, -t69, 0, 0; t69, t70, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:06
	% EndTime: 2020-11-04 21:34:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t74 = cos(qJ(1));
	t73 = sin(qJ(1));
	t72 = cos(pkin(9));
	t71 = sin(pkin(9));
	t1 = [t74 * t72, -t74 * t71, t73, t74 * pkin(1) + t73 * qJ(2) + 0; t73 * t72, -t73 * t71, -t74, t73 * pkin(1) - t74 * qJ(2) + 0; t71, t72, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:06
	% EndTime: 2020-11-04 21:34:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t81 = cos(qJ(1));
	t80 = sin(qJ(1));
	t79 = pkin(7) + qJ(2);
	t78 = pkin(9) + qJ(3);
	t77 = cos(t78);
	t76 = sin(t78);
	t75 = cos(pkin(9)) * pkin(2) + pkin(1);
	t1 = [t81 * t77, -t81 * t76, t80, t81 * t75 + t79 * t80 + 0; t80 * t77, -t80 * t76, -t81, t80 * t75 - t81 * t79 + 0; t76, t77, 0, sin(pkin(9)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:06
	% EndTime: 2020-11-04 21:34:06
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t86 = sin(pkin(10));
	t89 = sin(qJ(1));
	t95 = t89 * t86;
	t87 = cos(pkin(10));
	t94 = t89 * t87;
	t90 = cos(qJ(1));
	t93 = t90 * t86;
	t92 = t90 * t87;
	t85 = pkin(9) + qJ(3);
	t83 = sin(t85);
	t84 = cos(t85);
	t91 = pkin(3) * t84 + qJ(4) * t83 + cos(pkin(9)) * pkin(2) + pkin(1);
	t88 = pkin(7) + qJ(2);
	t1 = [t84 * t92 + t95, -t84 * t93 + t94, t90 * t83, t88 * t89 + t91 * t90 + 0; t84 * t94 - t93, -t84 * t95 - t92, t89 * t83, -t90 * t88 + t91 * t89 + 0; t83 * t87, -t83 * t86, -t84, t83 * pkin(3) - t84 * qJ(4) + sin(pkin(9)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:06
	% EndTime: 2020-11-04 21:34:06
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (52->24), mult. (53->28), div. (0->0), fcn. (70->8), ass. (0->18)
	t104 = sin(pkin(10));
	t107 = sin(qJ(1));
	t113 = t107 * t104;
	t105 = cos(pkin(10));
	t112 = t107 * t105;
	t108 = cos(qJ(1));
	t111 = t108 * t104;
	t110 = t108 * t105;
	t103 = pkin(9) + qJ(3);
	t101 = sin(t103);
	t102 = cos(t103);
	t109 = pkin(3) * t102 + qJ(4) * t101 + cos(pkin(9)) * pkin(2) + pkin(1);
	t106 = pkin(7) + qJ(2);
	t99 = t102 * t110 + t113;
	t98 = t102 * t111 - t112;
	t97 = t102 * t112 - t111;
	t96 = t102 * t113 + t110;
	t1 = [t99, t108 * t101, t98, t99 * pkin(4) + t98 * qJ(5) + t106 * t107 + t109 * t108 + 0; t97, t107 * t101, t96, t97 * pkin(4) + t96 * qJ(5) - t108 * t106 + t109 * t107 + 0; t101 * t105, -t102, t101 * t104, sin(pkin(9)) * pkin(2) - t102 * qJ(4) + pkin(6) + 0 + (pkin(4) * t105 + qJ(5) * t104 + pkin(3)) * t101; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:34:06
	% EndTime: 2020-11-04 21:34:06
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (84->35), mult. (80->40), div. (0->0), fcn. (96->16), ass. (0->25)
	t123 = sin(pkin(10));
	t125 = cos(pkin(10));
	t128 = sin(qJ(6));
	t130 = cos(qJ(6));
	t115 = t130 * t123 - t128 * t125;
	t129 = sin(qJ(1));
	t137 = t129 * t115;
	t116 = t128 * t123 + t130 * t125;
	t136 = t129 * t116;
	t131 = cos(qJ(1));
	t135 = t131 * t115;
	t134 = t131 * t116;
	t122 = pkin(9) + qJ(3);
	t132 = pkin(4) + pkin(5);
	t133 = qJ(5) * t125 - t132 * t123 - pkin(7) - qJ(2);
	t127 = qJ(4) - pkin(8);
	t126 = cos(pkin(9));
	t124 = sin(pkin(9));
	t121 = cos(t122);
	t120 = sin(t122);
	t119 = -pkin(10) + t122;
	t118 = pkin(10) + t122;
	t117 = qJ(5) * t123 + t132 * t125 + pkin(3);
	t114 = (t117 * t126 + t124 * t127) * cos(qJ(3)) + (-t124 * t117 + t127 * t126) * sin(qJ(3)) + t126 * pkin(2) + pkin(1);
	t1 = [t121 * t134 + t137, t121 * t135 - t136, -t131 * t120, t114 * t131 - t133 * t129 + 0; t121 * t136 - t135, t121 * t137 + t134, -t129 * t120, t114 * t129 + t133 * t131 + 0; t120 * t116, t120 * t115, t121, -t127 * t121 + t124 * pkin(2) + t120 * pkin(3) + 0 + pkin(6) + (sin(t119) / 0.2e1 + sin(t118) / 0.2e1) * t132 + (cos(t119) / 0.2e1 - cos(t118) / 0.2e1) * qJ(5); 0, 0, 0, 1;];
	Tc_mdh = t1;
end