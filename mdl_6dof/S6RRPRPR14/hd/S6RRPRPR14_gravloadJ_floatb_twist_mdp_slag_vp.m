% Calculate Gravitation load on the joints for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR14_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR14_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:36:57
% EndTime: 2019-03-09 11:37:01
% DurationCPUTime: 0.95s
% Computational Cost: add. (314->109), mult. (770->169), div. (0->0), fcn. (891->10), ass. (0->47)
t157 = MDP(12) - MDP(9) - MDP(22);
t156 = MDP(20) - MDP(23);
t154 = MDP(10) - MDP(13);
t151 = MDP(21) - MDP(24);
t114 = sin(qJ(2));
t115 = sin(qJ(1));
t118 = cos(qJ(2));
t119 = cos(qJ(1));
t145 = cos(pkin(6));
t134 = t119 * t145;
t100 = t114 * t134 + t115 * t118;
t135 = t115 * t145;
t102 = -t114 * t135 + t118 * t119;
t152 = -g(1) * t102 - g(2) * t100;
t112 = sin(qJ(6));
t116 = cos(qJ(6));
t113 = sin(qJ(4));
t117 = cos(qJ(4));
t111 = sin(pkin(6));
t139 = t111 * t119;
t99 = t114 * t115 - t118 * t134;
t124 = t113 * t139 + t99 * t117;
t150 = -t100 * t116 + t112 * t124;
t149 = t100 * t112 + t116 * t124;
t146 = g(3) * t111;
t142 = t111 * t114;
t141 = t111 * t115;
t140 = t111 * t118;
t138 = t112 * t117;
t137 = t116 * t117;
t136 = g(3) * (pkin(2) * t140 + qJ(3) * t142);
t128 = pkin(4) * t113 - qJ(5) * t117;
t101 = t119 * t114 + t118 * t135;
t127 = t119 * pkin(1) + t102 * pkin(2) + pkin(8) * t141 + qJ(3) * t101;
t125 = -t99 * t113 + t117 * t139;
t84 = -t101 * t117 + t113 * t141;
t97 = t145 * t113 + t117 * t140;
t123 = g(1) * t84 - g(2) * t124 + g(3) * t97;
t121 = -t115 * pkin(1) - t100 * pkin(2) + pkin(8) * t139 - t99 * qJ(3);
t81 = -g(1) * t101 - g(2) * t99 + g(3) * t140;
t98 = -t113 * t140 + t145 * t117;
t95 = t101 * pkin(2);
t93 = t99 * pkin(2);
t85 = t101 * t113 + t117 * t141;
t78 = t102 * t116 + t112 * t84;
t77 = -t102 * t112 + t116 * t84;
t1 = [(g(1) * t115 - g(2) * t119) * MDP(2) + (-g(1) * t121 - g(2) * t127) * MDP(14) + (-g(1) * (pkin(3) * t139 + pkin(4) * t125 - t100 * pkin(9) + qJ(5) * t124 + t121) - g(2) * (pkin(3) * t141 + pkin(4) * t85 + pkin(9) * t102 + qJ(5) * t84 + t127)) * MDP(25) + (-g(1) * t150 - g(2) * t78) * MDP(31) + (-g(1) * t149 - g(2) * t77) * MDP(32) + t151 * (g(1) * t124 + g(2) * t84) - t156 * (g(1) * t125 + g(2) * t85) - t154 * (g(1) * t99 - g(2) * t101) - t157 * (g(1) * t100 - g(2) * t102) + (-t111 * MDP(11) + MDP(3)) * (g(1) * t119 + g(2) * t115); (-g(1) * (qJ(3) * t102 - t95) - g(2) * (qJ(3) * t100 - t93) - t136) * MDP(14) + (-g(1) * (-pkin(9) * t101 - t95) - g(2) * (-pkin(9) * t99 - t93) - t136 - (pkin(9) * t118 + t128 * t114) * t146 + t152 * (qJ(3) + t128)) * MDP(25) + (-g(1) * (-t101 * t116 - t102 * t138) - g(2) * (-t100 * t138 - t116 * t99) - (-t114 * t138 + t116 * t118) * t146) * MDP(31) + (-g(1) * (t101 * t112 - t102 * t137) - g(2) * (-t100 * t137 + t112 * t99) - (-t112 * t118 - t114 * t137) * t146) * MDP(32) + t157 * t81 + (-t113 * t156 - t117 * t151 + t154) * (g(3) * t142 - t152); (MDP(14) + MDP(25)) * t81; (-g(1) * (-pkin(4) * t84 + qJ(5) * t85) - g(2) * (pkin(4) * t124 - qJ(5) * t125) - g(3) * (-pkin(4) * t97 + qJ(5) * t98)) * MDP(25) + (-MDP(31) * t112 - MDP(32) * t116 + t151) * (g(1) * t85 - g(2) * t125 + g(3) * t98) + t156 * t123; -t123 * MDP(25); (-g(1) * t77 + g(2) * t149 - g(3) * (-t112 * t142 + t116 * t97)) * MDP(31) + (g(1) * t78 - g(2) * t150 - g(3) * (-t112 * t97 - t116 * t142)) * MDP(32);];
taug  = t1;
