% Calculate Gravitation load on the joints for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:10:58
% EndTime: 2019-03-08 21:11:00
% DurationCPUTime: 0.66s
% Computational Cost: add. (268->88), mult. (677->142), div. (0->0), fcn. (784->10), ass. (0->48)
t143 = -MDP(10) - MDP(12) + MDP(17);
t107 = sin(qJ(3));
t142 = qJ(4) * t107 + pkin(2);
t140 = MDP(11) - MDP(14) - MDP(16);
t104 = sin(pkin(6));
t139 = g(3) * t104;
t138 = pkin(8) - qJ(5);
t110 = cos(qJ(3));
t108 = sin(qJ(2));
t111 = cos(qJ(2));
t105 = cos(pkin(10));
t135 = cos(pkin(6));
t117 = t105 * t135;
t134 = sin(pkin(10));
t90 = -t108 * t134 + t111 * t117;
t137 = t110 * t90;
t115 = t135 * t134;
t92 = -t105 * t108 - t111 * t115;
t136 = t110 * t92;
t132 = t104 * t108;
t131 = t104 * t110;
t130 = t104 * t111;
t106 = sin(qJ(6));
t129 = t106 * t107;
t109 = cos(qJ(6));
t128 = t107 * t109;
t127 = t107 * t111;
t126 = t109 * t111;
t125 = t110 * t111;
t124 = MDP(15) + MDP(19);
t123 = pkin(3) * t137 + t142 * t90;
t122 = pkin(3) * t136 + t142 * t92;
t91 = t108 * t117 + t111 * t134;
t79 = t105 * t131 + t107 * t91;
t80 = -t104 * t105 * t107 + t110 * t91;
t121 = -pkin(3) * t79 + qJ(4) * t80;
t118 = t104 * t134;
t93 = t105 * t111 - t108 * t115;
t81 = t107 * t93 - t110 * t118;
t82 = t107 * t118 + t110 * t93;
t120 = -pkin(3) * t81 + qJ(4) * t82;
t94 = t107 * t132 - t110 * t135;
t95 = t107 * t135 + t108 * t131;
t119 = -pkin(3) * t94 + qJ(4) * t95;
t116 = pkin(2) * t130 + pkin(8) * t132 + (pkin(3) * t125 + qJ(4) * t127) * t104;
t114 = g(1) * t81 + g(2) * t79 + g(3) * t94;
t112 = g(1) * t92 + g(2) * t90 + g(3) * t130;
t1 = [(-MDP(1) - t124) * g(3); (-g(1) * (pkin(8) * t93 + t122) - g(2) * (pkin(8) * t91 + t123) - g(3) * t116) * MDP(15) + (-g(1) * (pkin(4) * t136 + t138 * t93 + t122) - g(2) * (pkin(4) * t137 + t138 * t91 + t123) - g(3) * ((pkin(4) * t125 - qJ(5) * t108) * t104 + t116)) * MDP(19) + (-g(1) * (-t106 * t93 + t128 * t92) - g(2) * (-t106 * t91 + t128 * t90) - (-t106 * t108 + t107 * t126) * t139) * MDP(25) + (-g(1) * (-t109 * t93 - t129 * t92) - g(2) * (-t109 * t91 - t129 * t90) - (-t106 * t127 - t108 * t109) * t139) * MDP(26) + (MDP(4) - MDP(13) + MDP(18)) * (g(1) * t93 + g(2) * t91 + g(3) * t132) + (t140 * t107 + t110 * t143 - MDP(3)) * t112; (-g(1) * t120 - g(2) * t121 - g(3) * t119) * MDP(15) + (-g(1) * (-pkin(4) * t81 + t120) - g(2) * (-pkin(4) * t79 + t121) - g(3) * (-pkin(4) * t94 + t119)) * MDP(19) + (-MDP(25) * t109 + MDP(26) * t106 + t140) * (g(1) * t82 + g(2) * t80 + g(3) * t95) - t143 * t114; -t124 * t114; -t112 * MDP(19); (-g(1) * (-t106 * t81 + t109 * t92) - g(2) * (-t106 * t79 + t109 * t90) - g(3) * (t104 * t126 - t106 * t94)) * MDP(25) + (-g(1) * (-t106 * t92 - t109 * t81) - g(2) * (-t106 * t90 - t109 * t79) - g(3) * (-t106 * t130 - t109 * t94)) * MDP(26);];
taug  = t1;
