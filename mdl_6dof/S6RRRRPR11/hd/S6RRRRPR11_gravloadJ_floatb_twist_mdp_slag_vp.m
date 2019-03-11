% Calculate Gravitation load on the joints for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:23:35
% EndTime: 2019-03-09 23:23:38
% DurationCPUTime: 0.80s
% Computational Cost: add. (469->125), mult. (909->205), div. (0->0), fcn. (1087->12), ass. (0->54)
t103 = sin(qJ(3));
t107 = cos(qJ(3));
t100 = sin(pkin(6));
t135 = cos(qJ(1));
t117 = t100 * t135;
t104 = sin(qJ(2));
t105 = sin(qJ(1));
t108 = cos(qJ(2));
t126 = cos(pkin(6));
t114 = t126 * t135;
t88 = t104 * t114 + t105 * t108;
t80 = -t103 * t117 + t88 * t107;
t87 = t105 * t104 - t108 * t114;
t99 = qJ(4) + pkin(12) + qJ(6);
t96 = sin(t99);
t97 = cos(t99);
t144 = t80 * t96 - t87 * t97;
t143 = t80 * t97 + t87 * t96;
t102 = sin(qJ(4));
t106 = cos(qJ(4));
t142 = t80 * t102 - t87 * t106;
t141 = t87 * t102 + t80 * t106;
t139 = MDP(17) - MDP(25);
t137 = g(2) * t87;
t136 = g(2) * t88;
t134 = g(3) * t100;
t122 = t100 * t108;
t124 = t100 * t105;
t116 = t105 * t126;
t90 = -t104 * t116 + t135 * t108;
t84 = t103 * t124 + t90 * t107;
t89 = t135 * t104 + t108 * t116;
t74 = -t84 * t96 + t89 * t97;
t75 = t84 * t97 + t89 * t96;
t123 = t100 * t107;
t86 = t126 * t103 + t104 * t123;
t131 = (-g(1) * t74 + g(2) * t144 - g(3) * (-t97 * t122 - t86 * t96)) * MDP(32) + (g(1) * t75 + g(2) * t143 - g(3) * (t96 * t122 - t86 * t97)) * MDP(33);
t130 = t107 * t96;
t129 = t107 * t97;
t125 = t100 * t104;
t121 = t102 * t107;
t120 = t106 * t107;
t119 = t107 * t108;
t118 = pkin(4) * t102 + pkin(9);
t76 = -t84 * t102 + t89 * t106;
t101 = -qJ(5) - pkin(10);
t98 = t106 * pkin(4) + pkin(3);
t113 = t101 * t103 - t107 * t98 - pkin(2);
t79 = t88 * t103 + t107 * t117;
t83 = t90 * t103 - t105 * t123;
t85 = t103 * t125 - t126 * t107;
t112 = g(1) * t83 + g(2) * t79 + g(3) * t85;
t77 = t89 * t102 + t84 * t106;
t1 = [(g(1) * t105 - g(2) * t135) * MDP(2) + (g(1) * t135 + g(2) * t105) * MDP(3) + (g(1) * t88 - g(2) * t90) * MDP(9) + (-g(1) * t87 + g(2) * t89) * MDP(10) + (g(1) * t80 - g(2) * t84) * MDP(16) + (g(1) * t141 - g(2) * t77) * MDP(23) + (-g(1) * t142 - g(2) * t76) * MDP(24) + (-g(1) * (-t105 * pkin(1) - t88 * pkin(2) + pkin(8) * t117 + t101 * t79 - t118 * t87 - t80 * t98) - g(2) * (t135 * pkin(1) + t90 * pkin(2) + pkin(8) * t124 - t83 * t101 + t118 * t89 + t84 * t98)) * MDP(26) + (g(1) * t143 - g(2) * t75) * MDP(32) + (-g(1) * t144 - g(2) * t74) * MDP(33) + t139 * (-g(1) * t79 + g(2) * t83); (g(1) * t90 + g(3) * t125 + t136) * MDP(10) + (-g(1) * (t90 * t102 - t89 * t120) - g(2) * (t88 * t102 - t87 * t120) - (t102 * t104 + t106 * t119) * t134) * MDP(23) + (-g(1) * (t90 * t106 + t89 * t121) - g(2) * (t88 * t106 + t87 * t121) - (-t102 * t119 + t104 * t106) * t134) * MDP(24) + (-g(1) * (t113 * t89 + t118 * t90) - t118 * t136 - t113 * t137 - (t118 * t104 - t113 * t108) * t134) * MDP(26) + (-g(1) * (-t89 * t129 + t90 * t96) - g(2) * (-t87 * t129 + t88 * t96) - (t104 * t96 + t97 * t119) * t134) * MDP(32) + (-g(1) * (t89 * t130 + t90 * t97) - g(2) * (t87 * t130 + t88 * t97) - (t104 * t97 - t96 * t119) * t134) * MDP(33) + (-t107 * MDP(16) + t139 * t103 - MDP(9)) * (-g(1) * t89 + g(3) * t122 - t137); (-g(1) * (-t84 * t101 - t83 * t98) - g(2) * (-t80 * t101 - t79 * t98) - g(3) * (-t86 * t101 - t85 * t98)) * MDP(26) + t139 * (g(1) * t84 + g(2) * t80 + g(3) * t86) + (t106 * MDP(23) - t102 * MDP(24) + t97 * MDP(32) - t96 * MDP(33) + MDP(16)) * t112; (g(1) * t77 + g(2) * t141 - g(3) * (t102 * t122 - t86 * t106)) * MDP(24) + t131 + (pkin(4) * MDP(26) + MDP(23)) * (g(2) * t142 - g(3) * (-t86 * t102 - t106 * t122) - g(1) * t76); -t112 * MDP(26); t131;];
taug  = t1;
