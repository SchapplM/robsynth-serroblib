% Calculate Gravitation load on the joints for
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRRPP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:36:59
% EndTime: 2019-03-09 04:37:01
% DurationCPUTime: 0.58s
% Computational Cost: add. (373->87), mult. (474->118), div. (0->0), fcn. (465->8), ass. (0->43)
t143 = MDP(18) - MDP(21) - MDP(24);
t141 = MDP(17) - MDP(20) + MDP(25);
t140 = MDP(11) - MDP(19) - MDP(23);
t105 = sin(qJ(3));
t103 = qJ(1) + pkin(9);
t98 = sin(t103);
t99 = cos(t103);
t119 = g(1) * t99 + g(2) * t98;
t139 = t119 * t105;
t136 = g(1) * t98;
t100 = t105 * pkin(8);
t108 = cos(qJ(3));
t101 = t108 * pkin(3);
t133 = -pkin(4) - qJ(6);
t132 = t105 * t99;
t131 = t108 * t98;
t130 = t108 * t99;
t104 = sin(qJ(4));
t129 = t104 * t105;
t128 = t104 * t108;
t107 = cos(qJ(4));
t127 = t105 * t107;
t126 = t107 * t108;
t125 = MDP(22) + MDP(26);
t124 = -pkin(2) - t101;
t123 = -qJ(5) * t104 - pkin(3);
t81 = t99 * t107 + t98 * t128;
t82 = -t99 * t104 + t98 * t126;
t122 = -t81 * pkin(4) + t82 * qJ(5);
t83 = -t98 * t107 + t99 * t128;
t84 = t98 * t104 + t99 * t126;
t121 = -t83 * pkin(4) + t84 * qJ(5);
t120 = pkin(4) * t126 + qJ(5) * t128 + t100 + t101;
t106 = sin(qJ(1));
t115 = -t106 * pkin(1) - t82 * pkin(4) + t99 * pkin(7) - t81 * qJ(5);
t113 = g(1) * t83 + g(2) * t81 + g(3) * t129;
t112 = g(1) * t84 + g(2) * t82 + g(3) * t127;
t109 = cos(qJ(1));
t111 = t109 * pkin(1) + t99 * pkin(2) + pkin(3) * t130 + t84 * pkin(4) + t98 * pkin(7) + pkin(8) * t132 + t83 * qJ(5);
t92 = qJ(5) * t127;
t88 = pkin(8) * t130;
t86 = pkin(8) * t131;
t1 = [(g(1) * t109 + g(2) * t106) * MDP(3) + (-g(1) * t115 - g(2) * t111 - (t124 - t100) * t136) * MDP(22) + (-g(1) * (-t82 * qJ(6) + t115) - g(2) * (pkin(5) * t132 + t84 * qJ(6) + t111) - ((-pkin(5) - pkin(8)) * t105 + t124) * t136) * MDP(26) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t106 - g(2) * t109) + t141 * (g(1) * t82 - g(2) * t84) - t143 * (g(1) * t81 - g(2) * t83) + (t108 * MDP(10) - t140 * t105) * (-g(2) * t99 + t136); (-MDP(4) - t125) * g(3); (-g(1) * t88 - g(2) * t86 - g(3) * t120 + (pkin(4) * t107 - t123) * t139) * MDP(22) + (-g(1) * (pkin(5) * t130 + t88) - g(2) * (pkin(5) * t131 + t86) - g(3) * (qJ(6) * t126 + t120) + (-g(3) * pkin(5) + t119 * (-t133 * t107 - t123)) * t105) * MDP(26) + t140 * (g(3) * t105 + t119 * t108) + (-t104 * t143 + t141 * t107 + MDP(10)) * (-g(3) * t108 + t139); (-g(1) * t121 - g(2) * t122 - g(3) * (-pkin(4) * t129 + t92)) * MDP(22) + (-g(1) * (-t83 * qJ(6) + t121) - g(2) * (-t81 * qJ(6) + t122) - g(3) * (t133 * t129 + t92)) * MDP(26) + t141 * t113 + t143 * t112; -t125 * t113; -t112 * MDP(26);];
taug  = t1;
