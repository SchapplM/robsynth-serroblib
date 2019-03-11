% Calculate Gravitation load on the joints for
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:46:31
% EndTime: 2019-03-09 19:46:34
% DurationCPUTime: 1.02s
% Computational Cost: add. (525->142), mult. (988->235), div. (0->0), fcn. (1193->14), ass. (0->54)
t139 = MDP(17) - MDP(20);
t105 = sin(qJ(3));
t107 = cos(qJ(3));
t103 = sin(pkin(6));
t136 = cos(qJ(1));
t118 = t103 * t136;
t106 = sin(qJ(2));
t108 = cos(qJ(2));
t124 = cos(pkin(6));
t114 = t124 * t136;
t135 = sin(qJ(1));
t88 = t106 * t114 + t135 * t108;
t80 = -t105 * t118 + t107 * t88;
t87 = t135 * t106 - t108 * t114;
t101 = pkin(12) + qJ(5);
t98 = sin(t101);
t99 = cos(t101);
t144 = t80 * t98 - t87 * t99;
t143 = t80 * t99 + t87 * t98;
t100 = qJ(6) + t101;
t96 = sin(t100);
t97 = cos(t100);
t142 = t80 * t96 - t87 * t97;
t141 = t80 * t97 + t87 * t96;
t113 = t124 * t135;
t89 = t136 * t106 + t108 * t113;
t140 = -g(1) * t89 - g(2) * t87;
t134 = g(3) * t103;
t121 = t103 * t108;
t117 = t103 * t135;
t90 = -t106 * t113 + t136 * t108;
t84 = t105 * t117 + t90 * t107;
t74 = -t84 * t96 + t89 * t97;
t75 = t84 * t97 + t89 * t96;
t122 = t103 * t106;
t86 = t124 * t105 + t107 * t122;
t129 = (-g(1) * t74 + g(2) * t142 - g(3) * (-t97 * t121 - t86 * t96)) * MDP(34) + (g(1) * t75 + g(2) * t141 - g(3) * (t96 * t121 - t86 * t97)) * MDP(35);
t128 = t107 * t96;
t127 = t107 * t97;
t126 = t107 * t98;
t125 = t107 * t99;
t102 = sin(pkin(12));
t123 = t102 * t107;
t104 = cos(pkin(12));
t120 = t104 * t107;
t119 = t107 * t108;
t115 = -g(1) * t90 - g(2) * t88;
t79 = t88 * t105 + t107 * t118;
t83 = t105 * t90 - t107 * t117;
t85 = t105 * t122 - t124 * t107;
t111 = g(1) * t83 + g(2) * t79 + g(3) * t85;
t77 = t84 * t99 + t89 * t98;
t76 = -t84 * t98 + t89 * t99;
t1 = [(g(1) * t135 - g(2) * t136) * MDP(2) + (g(1) * t136 + g(2) * t135) * MDP(3) + (g(1) * t88 - g(2) * t90) * MDP(9) + (-g(1) * t87 + g(2) * t89) * MDP(10) + (g(1) * t80 - g(2) * t84) * MDP(16) + (-g(1) * (-t102 * t87 - t104 * t80) - g(2) * (t102 * t89 + t104 * t84)) * MDP(18) + (-g(1) * (t102 * t80 - t104 * t87) - g(2) * (-t102 * t84 + t104 * t89)) * MDP(19) + (-g(1) * (-t135 * pkin(1) - t88 * pkin(2) - pkin(3) * t80 + pkin(8) * t118 - t87 * pkin(9) - qJ(4) * t79) - g(2) * (t136 * pkin(1) + t90 * pkin(2) + t84 * pkin(3) + pkin(8) * t117 + t89 * pkin(9) + t83 * qJ(4))) * MDP(21) + (g(1) * t143 - g(2) * t77) * MDP(27) + (-g(1) * t144 - g(2) * t76) * MDP(28) + (g(1) * t141 - g(2) * t75) * MDP(34) + (-g(1) * t142 - g(2) * t74) * MDP(35) + t139 * (-g(1) * t79 + g(2) * t83); (g(3) * t122 - t115) * MDP(10) + (-g(1) * (t102 * t90 - t89 * t120) - g(2) * (t102 * t88 - t87 * t120) - (t102 * t106 + t104 * t119) * t134) * MDP(18) + (-g(1) * (t104 * t90 + t89 * t123) - g(2) * (t104 * t88 + t87 * t123) - (-t102 * t119 + t104 * t106) * t134) * MDP(19) + ((-t106 * t134 + t115) * pkin(9) + (-t108 * t134 - t140) * (pkin(3) * t107 + qJ(4) * t105 + pkin(2))) * MDP(21) + (-g(1) * (-t89 * t125 + t90 * t98) - g(2) * (-t87 * t125 + t88 * t98) - (t106 * t98 + t99 * t119) * t134) * MDP(27) + (-g(1) * (t89 * t126 + t90 * t99) - g(2) * (t87 * t126 + t88 * t99) - (t106 * t99 - t98 * t119) * t134) * MDP(28) + (-g(1) * (-t89 * t127 + t90 * t96) - g(2) * (-t87 * t127 + t88 * t96) - (t106 * t96 + t97 * t119) * t134) * MDP(34) + (-g(1) * (t89 * t128 + t90 * t97) - g(2) * (t87 * t128 + t88 * t97) - (t106 * t97 - t96 * t119) * t134) * MDP(35) + (-t107 * MDP(16) + t139 * t105 - MDP(9)) * (g(3) * t121 + t140); (-g(1) * (-pkin(3) * t83 + qJ(4) * t84) - g(2) * (-pkin(3) * t79 + qJ(4) * t80) - g(3) * (-pkin(3) * t85 + qJ(4) * t86)) * MDP(21) + t139 * (g(1) * t84 + g(2) * t80 + g(3) * t86) + (MDP(18) * t104 - MDP(19) * t102 + MDP(27) * t99 - MDP(28) * t98 + MDP(34) * t97 - MDP(35) * t96 + MDP(16)) * t111; -t111 * MDP(21); (-g(1) * t76 + g(2) * t144 - g(3) * (-t99 * t121 - t86 * t98)) * MDP(27) + (g(1) * t77 + g(2) * t143 - g(3) * (t98 * t121 - t86 * t99)) * MDP(28) + t129; t129;];
taug  = t1;
