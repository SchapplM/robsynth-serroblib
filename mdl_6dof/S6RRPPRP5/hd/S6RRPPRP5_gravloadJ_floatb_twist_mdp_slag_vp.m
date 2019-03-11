% Calculate Gravitation load on the joints for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:43:58
% EndTime: 2019-03-09 08:44:00
% DurationCPUTime: 0.72s
% Computational Cost: add. (281->98), mult. (449->127), div. (0->0), fcn. (417->8), ass. (0->45)
t100 = sin(qJ(1));
t102 = cos(qJ(1));
t129 = -g(1) * t102 - g(2) * t100;
t99 = sin(qJ(2));
t132 = t129 * t99;
t131 = MDP(10) - MDP(13);
t128 = MDP(24) + MDP(26);
t127 = MDP(25) - MDP(28);
t130 = MDP(9) - MDP(12) + MDP(17) + MDP(27);
t101 = cos(qJ(2));
t89 = t99 * qJ(3);
t119 = t101 * pkin(2) + t89;
t125 = g(3) * t99;
t77 = -t129 * t101 + t125;
t96 = sin(pkin(9));
t126 = pkin(4) * t96;
t124 = g(1) * t100;
t121 = g(3) * t101;
t118 = t100 * t99;
t98 = -pkin(8) - qJ(4);
t117 = t101 * t98;
t116 = t102 * t99;
t115 = qJ(4) * t101;
t114 = t101 * t102;
t113 = -MDP(18) - MDP(29);
t112 = t96 * t116;
t111 = pkin(2) * t114 + t100 * pkin(7) + (pkin(1) + t89) * t102;
t82 = t100 * t101 * qJ(3);
t84 = qJ(3) * t114;
t109 = -g(1) * t84 - g(2) * t82;
t107 = -pkin(1) - t119;
t95 = pkin(9) + qJ(5);
t87 = sin(t95);
t88 = cos(t95);
t106 = pkin(5) * t87 - qJ(6) * t88 + t126;
t72 = t100 * t87 - t88 * t116;
t74 = t102 * t87 + t88 * t118;
t68 = g(1) * t72 - g(2) * t74 + t88 * t121;
t97 = cos(pkin(9));
t92 = t102 * pkin(7);
t86 = pkin(4) * t97 + pkin(3);
t76 = -t121 - t132;
t75 = t102 * t88 - t87 * t118;
t73 = t100 * t88 + t87 * t116;
t1 = [(-g(1) * t92 - g(2) * t111 - t107 * t124) * MDP(14) + (-g(1) * (t102 * t97 - t96 * t118) - g(2) * (t100 * t97 + t112)) * MDP(15) + (-g(1) * (-t102 * t96 - t97 * t118) - g(2) * (-t100 * t96 + t97 * t116)) * MDP(16) + (-g(1) * (pkin(3) * t102 + t92) - g(2) * (qJ(4) * t114 + t111) + (-g(1) * (t107 - t115) - g(2) * pkin(3)) * t100) * MDP(18) + (-g(1) * (t75 * pkin(5) + t74 * qJ(6) + t102 * t86 + t92) - g(2) * (pkin(4) * t112 + t73 * pkin(5) + t72 * qJ(6) - t98 * t114 + t111) + (-g(1) * (-t99 * t126 + t107 + t117) - g(2) * t86) * t100) * MDP(29) - (MDP(3) - MDP(11)) * t129 + t128 * (-g(1) * t75 - g(2) * t73) + t127 * (g(1) * t74 + g(2) * t72) + (t130 * t101 - t131 * t99 + MDP(2)) * (-g(2) * t102 + t124); (-g(1) * (-pkin(2) * t116 + t84) - g(2) * (-pkin(2) * t118 + t82) - g(3) * t119) * MDP(14) + (-g(3) * (t115 + t119) + t109 - (pkin(2) + qJ(4)) * t132) * MDP(18) + (-g(3) * (-t117 + t119) - t106 * t125 + t109 + t129 * ((-pkin(2) + t98) * t99 + t106 * t101)) * MDP(29) + t130 * t76 + (-t96 * MDP(15) - t97 * MDP(16) - t127 * t88 - t128 * t87 + t131) * t77; (-MDP(14) + t113) * t76; t113 * t77; (-g(1) * (-pkin(5) * t72 + qJ(6) * t73) - g(2) * (pkin(5) * t74 - qJ(6) * t75) - (-pkin(5) * t88 - qJ(6) * t87) * t121) * MDP(29) + t128 * t68 - t127 * (-g(1) * t73 + g(2) * t75 + t87 * t121); -t68 * MDP(29);];
taug  = t1;
