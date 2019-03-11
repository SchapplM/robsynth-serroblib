% Calculate Gravitation load on the joints for
% S6RPRRPP2
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
%   see S6RPRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRRPP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:33:22
% EndTime: 2019-03-09 04:33:23
% DurationCPUTime: 0.48s
% Computational Cost: add. (369->87), mult. (469->120), div. (0->0), fcn. (458->8), ass. (0->43)
t140 = MDP(11) - MDP(20) + MDP(25);
t138 = MDP(17) + MDP(19) + MDP(23);
t137 = MDP(18) - MDP(21) - MDP(24);
t107 = sin(qJ(3));
t105 = qJ(1) + pkin(9);
t100 = sin(t105);
t101 = cos(t105);
t119 = g(1) * t101 + g(2) * t100;
t139 = t119 * t107;
t136 = -pkin(4) - pkin(5);
t135 = g(1) * t100;
t102 = t107 * pkin(8);
t110 = cos(qJ(3));
t103 = t110 * pkin(3);
t132 = t100 * t110;
t131 = t101 * t107;
t130 = t101 * t110;
t106 = sin(qJ(4));
t129 = t106 * t107;
t128 = t106 * t110;
t109 = cos(qJ(4));
t127 = t107 * t109;
t126 = t109 * t110;
t125 = MDP(22) + MDP(26);
t124 = -pkin(2) - t103;
t123 = -qJ(5) * t106 - pkin(3);
t82 = t100 * t128 + t101 * t109;
t83 = t100 * t126 - t101 * t106;
t122 = -t82 * pkin(4) + qJ(5) * t83;
t84 = -t100 * t109 + t101 * t128;
t85 = t100 * t106 + t101 * t126;
t121 = -t84 * pkin(4) + qJ(5) * t85;
t120 = pkin(4) * t126 + qJ(5) * t128 + t102 + t103;
t108 = sin(qJ(1));
t116 = -pkin(1) * t108 - t83 * pkin(4) + t101 * pkin(7) - qJ(5) * t82;
t114 = g(1) * t84 + g(2) * t82 + g(3) * t129;
t111 = cos(qJ(1));
t112 = t111 * pkin(1) + t101 * pkin(2) + pkin(3) * t130 + t85 * pkin(4) + t100 * pkin(7) + pkin(8) * t131 + qJ(5) * t84;
t76 = -g(3) * t110 + t139;
t94 = qJ(5) * t127;
t89 = pkin(8) * t130;
t87 = pkin(8) * t132;
t1 = [(g(1) * t111 + g(2) * t108) * MDP(3) + (-g(2) * t101 + t135) * t110 * MDP(10) + (-g(1) * t116 - g(2) * t112 - (t124 - t102) * t135) * MDP(22) + (-g(1) * (-pkin(5) * t83 + t116) - g(2) * (pkin(5) * t85 - qJ(6) * t131 + t112) - ((-pkin(8) + qJ(6)) * t107 + t124) * t135) * MDP(26) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t108 - g(2) * t111) - t140 * (-g(2) * t131 + t107 * t135) + t138 * (g(1) * t83 - g(2) * t85) - t137 * (g(1) * t82 - g(2) * t84); (-MDP(4) - t125) * g(3); (-g(1) * t89 - g(2) * t87 - g(3) * t120 + (pkin(4) * t109 - t123) * t139) * MDP(22) + (-g(1) * (-qJ(6) * t130 + t89) - g(2) * (-qJ(6) * t132 + t87) - g(3) * (pkin(5) * t126 + t120) + (g(3) * qJ(6) + t119 * (-t136 * t109 - t123)) * t107) * MDP(26) + t140 * (g(3) * t107 + t119 * t110) + (-t137 * t106 + t138 * t109 + MDP(10)) * t76; (-g(1) * t121 - g(2) * t122 - g(3) * (-pkin(4) * t129 + t94)) * MDP(22) + (-g(1) * (-pkin(5) * t84 + t121) - g(2) * (-pkin(5) * t82 + t122) - g(3) * (t136 * t129 + t94)) * MDP(26) + t138 * t114 + t137 * (g(1) * t85 + g(2) * t83 + g(3) * t127); -t125 * t114; t76 * MDP(26);];
taug  = t1;
