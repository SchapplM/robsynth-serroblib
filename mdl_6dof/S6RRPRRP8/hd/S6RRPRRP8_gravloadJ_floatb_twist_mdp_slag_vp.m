% Calculate Gravitation load on the joints for
% S6RRPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRRP8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:25:16
% EndTime: 2019-03-09 12:25:18
% DurationCPUTime: 0.63s
% Computational Cost: add. (476->97), mult. (493->139), div. (0->0), fcn. (482->10), ass. (0->44)
t142 = MDP(27) + MDP(29);
t143 = MDP(28) - MDP(31);
t112 = sin(qJ(1));
t114 = cos(qJ(1));
t122 = g(1) * t114 + g(2) * t112;
t108 = pkin(10) + qJ(4);
t100 = sin(t108);
t111 = sin(qJ(2));
t133 = g(3) * t111;
t101 = cos(t108);
t113 = cos(qJ(2));
t128 = t112 * t113;
t86 = t100 * t128 + t101 * t114;
t127 = t113 * t114;
t88 = -t100 * t127 + t112 * t101;
t140 = -g(1) * t88 + g(2) * t86 + t100 * t133;
t139 = MDP(10) - MDP(13) - MDP(30);
t90 = -g(3) * t113 + t111 * t122;
t136 = g(1) * t112;
t102 = qJ(5) + t108;
t98 = sin(t102);
t131 = t111 * t98;
t99 = cos(t102);
t130 = t111 * t99;
t107 = -pkin(9) - pkin(8) - qJ(3);
t129 = t107 * t111;
t126 = t114 * pkin(1) + t112 * pkin(7);
t82 = t114 * t99 + t128 * t98;
t84 = -t112 * t99 + t127 * t98;
t76 = g(1) * t84 + g(2) * t82 + g(3) * t131;
t83 = -t114 * t98 + t128 * t99;
t85 = t112 * t98 + t127 * t99;
t124 = t142 * t76 + t143 * (g(1) * t85 + g(2) * t83 + g(3) * t130);
t120 = pkin(2) * t113 + qJ(3) * t111;
t110 = cos(pkin(10));
t93 = pkin(3) * t110 + pkin(4) * t101 + pkin(2);
t118 = pkin(5) * t99 + qJ(6) * t98 + t93;
t115 = -g(1) * (-t84 * pkin(5) + t85 * qJ(6)) - g(2) * (-t82 * pkin(5) + t83 * qJ(6)) - g(3) * (-pkin(5) * t131 + qJ(6) * t130);
t109 = sin(pkin(10));
t104 = t114 * pkin(7);
t94 = pkin(3) * t109 + pkin(4) * t100;
t89 = t112 * t100 + t101 * t127;
t87 = t100 * t114 - t101 * t128;
t1 = [t122 * MDP(3) + (-g(1) * (t109 * t114 - t110 * t128) - g(2) * (t112 * t109 + t110 * t127)) * MDP(11) + (-g(1) * (t109 * t128 + t110 * t114) - g(2) * (-t109 * t127 + t112 * t110)) * MDP(12) + (-g(1) * t104 - g(2) * (t114 * t120 + t126) - (-pkin(1) - t120) * t136) * MDP(14) + (-g(1) * t87 - g(2) * t89) * MDP(20) + (-g(1) * t86 - g(2) * t88) * MDP(21) + (-g(1) * (-t83 * pkin(5) - t82 * qJ(6) + t114 * t94 + t104) - g(2) * (t85 * pkin(5) + t84 * qJ(6) - t114 * t129 + t127 * t93 + t126) + (-g(1) * (-t113 * t93 - pkin(1) + t129) - g(2) * t94) * t112) * MDP(32) + t142 * (g(1) * t83 - g(2) * t85) - t143 * (g(1) * t82 - g(2) * t84) + (t113 * MDP(9) - t139 * t111 + MDP(2)) * (-g(2) * t114 + t136); (-g(3) * t120 + t122 * (pkin(2) * t111 - qJ(3) * t113)) * MDP(14) + ((-g(3) * t118 + t107 * t122) * t113 + (g(3) * t107 + t122 * t118) * t111) * MDP(32) + t139 * (t113 * t122 + t133) + (t110 * MDP(11) - t109 * MDP(12) + t101 * MDP(20) - t100 * MDP(21) + t142 * t99 - t143 * t98 + MDP(9)) * t90; (-MDP(14) - MDP(32)) * t90; t140 * MDP(20) + (g(1) * t89 - g(2) * t87 + t101 * t133) * MDP(21) + (t140 * pkin(4) + t115) * MDP(32) + t124; MDP(32) * t115 + t124; -t76 * MDP(32);];
taug  = t1;
