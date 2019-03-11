% Calculate Gravitation load on the joints for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:26:09
% EndTime: 2019-03-09 09:26:11
% DurationCPUTime: 0.62s
% Computational Cost: add. (259->89), mult. (535->135), div. (0->0), fcn. (578->10), ass. (0->50)
t109 = cos(qJ(2));
t106 = sin(qJ(2));
t107 = sin(qJ(1));
t110 = cos(qJ(1));
t125 = g(1) * t110 + g(2) * t107;
t140 = t125 * t106;
t81 = -g(3) * t109 + t140;
t143 = MDP(11) + MDP(15);
t142 = -MDP(12) + MDP(17);
t141 = MDP(10) - MDP(13) - MDP(16);
t132 = t109 * pkin(2) + t106 * qJ(3);
t138 = g(1) * t107;
t135 = g(3) * t106;
t103 = sin(pkin(10));
t104 = cos(pkin(10));
t102 = qJ(5) + qJ(6);
t94 = sin(t102);
t95 = cos(t102);
t120 = t103 * t95 - t104 * t94;
t113 = g(3) * t120;
t119 = t103 * t94 + t104 * t95;
t130 = t107 * t109;
t83 = t103 * t130 + t104 * t110;
t128 = t110 * t103;
t84 = t104 * t130 - t128;
t122 = t83 * t94 + t84 * t95;
t123 = t83 * t95 - t84 * t94;
t85 = -t107 * t104 + t109 * t128;
t129 = t109 * t110;
t86 = t107 * t103 + t104 * t129;
t74 = t85 * t95 - t86 * t94;
t75 = t85 * t94 + t86 * t95;
t133 = (-g(1) * t74 - g(2) * t123 - t106 * t113) * MDP(31) + (g(1) * t75 + g(2) * t122 + t119 * t135) * MDP(32);
t131 = t106 * t110;
t127 = t110 * pkin(1) + pkin(2) * t129 + t107 * pkin(7) + qJ(3) * t131;
t121 = pkin(3) * t104 + qJ(4) * t103;
t105 = sin(qJ(5));
t108 = cos(qJ(5));
t118 = t105 * t84 - t108 * t83;
t117 = t105 * t83 + t108 * t84;
t116 = t103 * t108 - t104 * t105;
t115 = t103 * t105 + t104 * t108;
t112 = g(3) * t116;
t111 = (-pkin(1) - t132) * t138;
t99 = t110 * pkin(7);
t90 = qJ(3) * t129;
t88 = qJ(3) * t130;
t77 = t105 * t85 + t108 * t86;
t76 = -t105 * t86 + t108 * t85;
t1 = [t125 * MDP(3) + (-g(1) * t99 - g(2) * t127 - t111) * MDP(14) + (-g(1) * (-pkin(3) * t84 - qJ(4) * t83 + t99) - g(2) * (pkin(3) * t86 + qJ(4) * t85 + t127) - t111) * MDP(18) + (g(1) * t117 - g(2) * t77) * MDP(24) + (-g(1) * t118 - g(2) * t76) * MDP(25) + (g(1) * t122 - g(2) * t75) * MDP(31) + (g(1) * t123 - g(2) * t74) * MDP(32) + t143 * (g(1) * t84 - g(2) * t86) + t142 * (g(1) * t83 - g(2) * t85) + (MDP(9) * t109 - t141 * t106 + MDP(2)) * (-g(2) * t110 + t138); (-g(1) * (-pkin(2) * t131 + t90) - g(2) * (-pkin(2) * t106 * t107 + t88) - g(3) * t132) * MDP(14) + (-g(1) * t90 - g(2) * t88 - g(3) * (t121 * t109 + t132) + (pkin(2) + t121) * t140) * MDP(18) + (-t109 * t112 + t116 * t140) * MDP(25) + (-t109 * t113 + t120 * t140) * MDP(32) + t141 * (t125 * t109 + t135) + (t115 * MDP(24) + t119 * MDP(31) + t142 * t103 + t143 * t104 + MDP(9)) * t81; (-MDP(14) - MDP(18)) * t81; (-g(1) * t85 - g(2) * t83 - t103 * t135) * MDP(18); (-g(1) * t76 + g(2) * t118 - t106 * t112) * MDP(24) + (g(1) * t77 + g(2) * t117 + t115 * t135) * MDP(25) + t133; t133;];
taug  = t1;
