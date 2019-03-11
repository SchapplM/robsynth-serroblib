% Calculate Gravitation load on the joints for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRRPP7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:52:12
% EndTime: 2019-03-09 04:52:14
% DurationCPUTime: 0.53s
% Computational Cost: add. (214->91), mult. (483->126), div. (0->0), fcn. (470->6), ass. (0->46)
t142 = MDP(20) - MDP(23) - MDP(26);
t144 = MDP(13) - MDP(22) + MDP(27);
t143 = MDP(19) + MDP(21) + MDP(25);
t141 = -pkin(1) - pkin(7);
t140 = -pkin(4) - pkin(5);
t109 = sin(qJ(1));
t139 = g(1) * t109;
t112 = cos(qJ(1));
t138 = g(2) * t112;
t108 = sin(qJ(3));
t137 = g(3) * t108;
t136 = -pkin(8) + qJ(6);
t107 = sin(qJ(4));
t111 = cos(qJ(3));
t135 = t107 * t111;
t134 = t108 * t109;
t133 = t108 * t112;
t132 = t109 * t107;
t110 = cos(qJ(4));
t131 = t109 * t110;
t130 = t109 * t111;
t129 = t110 * t111;
t128 = t111 * t112;
t127 = t112 * t110;
t126 = t112 * pkin(1) + t109 * qJ(2);
t125 = MDP(24) + MDP(28);
t124 = t107 * t130;
t123 = t109 * t129;
t122 = -qJ(5) * t107 - pkin(3);
t84 = t108 * t132 - t127;
t85 = t107 * t112 + t108 * t131;
t121 = -t84 * pkin(4) + qJ(5) * t85;
t86 = t107 * t133 + t131;
t87 = t108 * t127 - t132;
t120 = t86 * pkin(4) - qJ(5) * t87;
t119 = pkin(3) * t130 + pkin(4) * t123 + pkin(8) * t134 + qJ(5) * t124;
t92 = -t138 + t139;
t117 = -pkin(4) * t110 + t122;
t116 = pkin(3) * t134 + t85 * pkin(4) + t112 * pkin(7) + t84 * qJ(5) + t126;
t115 = g(1) * t84 - g(2) * t86 + g(3) * t135;
t102 = t112 * qJ(2);
t113 = pkin(3) * t133 + t87 * pkin(4) - pkin(8) * t128 + t86 * qJ(5) + t102;
t79 = -t92 * t111 + t137;
t103 = t111 * pkin(8);
t94 = qJ(5) * t129;
t1 = [(-g(1) * (-t109 * pkin(1) + t102) - g(2) * t126) * MDP(6) + (-g(1) * (t141 * t109 + t113) - g(2) * (-pkin(8) * t130 + t116)) * MDP(24) + (-g(1) * (t87 * pkin(5) + qJ(6) * t128 + t113) - g(2) * (pkin(5) * t85 + t116) + (-g(2) * t136 * t111 - g(1) * t141) * t109) * MDP(28) + (MDP(2) - MDP(4)) * t92 + t143 * (-g(1) * t87 - g(2) * t85) + (-MDP(12) * t108 - t144 * t111 + MDP(3) - MDP(5)) * (g(1) * t112 + g(2) * t109) + t142 * (g(1) * t86 + g(2) * t84); (-MDP(6) - t125) * t92; (-g(1) * t119 - g(3) * t103 - t117 * t137 - (-pkin(8) * t108 + t117 * t111) * t138) * MDP(24) + (-g(1) * (pkin(5) * t123 + t119) - g(3) * (-qJ(6) * t111 + t103) + (qJ(6) * t139 - g(3) * (-pkin(5) * t110 + t117)) * t108 - (t136 * t108 + (t140 * t110 + t122) * t111) * t138) * MDP(28) + t144 * (g(3) * t111 + t92 * t108) + t142 * (g(1) * t124 + (-g(2) * t128 - t137) * t107) + (t143 * t110 + MDP(12)) * t79; (-g(1) * t121 - g(2) * t120 - g(3) * (-pkin(4) * t135 + t94)) * MDP(24) + (-g(1) * (-pkin(5) * t84 + t121) - g(2) * (pkin(5) * t86 + t120) - g(3) * (t140 * t135 + t94)) * MDP(28) + t143 * t115 + t142 * (g(1) * t85 - g(2) * t87 + g(3) * t129); -t125 * t115; t79 * MDP(28);];
taug  = t1;
