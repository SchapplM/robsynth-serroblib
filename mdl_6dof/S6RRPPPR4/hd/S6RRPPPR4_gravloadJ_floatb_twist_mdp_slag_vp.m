% Calculate Gravitation load on the joints for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:19:45
% EndTime: 2019-03-09 08:19:47
% DurationCPUTime: 0.66s
% Computational Cost: add. (196->81), mult. (455->111), div. (0->0), fcn. (445->8), ass. (0->46)
t144 = MDP(18) + MDP(22);
t104 = sin(qJ(2));
t107 = cos(qJ(2));
t105 = sin(qJ(1));
t108 = cos(qJ(1));
t84 = g(1) * t108 + g(2) * t105;
t77 = g(3) * t104 + t107 * t84;
t143 = t104 * t84;
t141 = MDP(10) - MDP(13);
t140 = MDP(15) + MDP(19);
t139 = MDP(16) - MDP(21);
t138 = MDP(9) - MDP(12) + MDP(17) + MDP(20);
t134 = g(1) * t105;
t130 = g(3) * t107;
t96 = t107 * pkin(2);
t129 = pkin(2) + qJ(4);
t92 = t104 * qJ(3);
t128 = t96 + t92;
t97 = t108 * pkin(7);
t127 = t108 * pkin(3) + t97;
t126 = t104 * t105;
t125 = t104 * t108;
t102 = cos(pkin(9));
t124 = t105 * t102;
t123 = t107 * t108;
t121 = -pkin(1) - t92;
t120 = t108 * pkin(1) + pkin(2) * t123 + t105 * pkin(7) + qJ(3) * t125;
t103 = sin(qJ(6));
t106 = cos(qJ(6));
t101 = sin(pkin(9));
t80 = t101 * t108 + t104 * t124;
t81 = -t101 * t126 + t102 * t108;
t116 = t103 * t81 - t106 * t80;
t115 = t103 * t80 + t106 * t81;
t114 = t105 * pkin(3) + qJ(4) * t123 + t120;
t113 = t101 * t106 - t102 * t103;
t112 = t101 * t103 + t102 * t106;
t109 = (-t107 * t129 + t121) * t134;
t90 = qJ(3) * t123;
t87 = t105 * t107 * qJ(3);
t79 = t101 * t125 + t124;
t78 = t101 * t105 - t102 * t125;
t76 = -t130 + t143;
t72 = t103 * t78 + t106 * t79;
t71 = -t103 * t79 + t106 * t78;
t1 = [(-g(1) * t97 - g(2) * t120 - (t121 - t96) * t134) * MDP(14) + (-g(1) * t127 - g(2) * t114 - t109) * MDP(18) + (-g(1) * (pkin(4) * t81 + qJ(5) * t80 + t127) - g(2) * (pkin(4) * t79 + qJ(5) * t78 + t114) - t109) * MDP(22) + (-g(1) * t115 - g(2) * t72) * MDP(28) + (g(1) * t116 - g(2) * t71) * MDP(29) + (MDP(3) - MDP(11)) * t84 + t140 * (-g(1) * t81 - g(2) * t79) + t139 * (g(1) * t80 + g(2) * t78) + (-t141 * t104 + t138 * t107 + MDP(2)) * (-g(2) * t108 + t134); (-g(1) * (-pkin(2) * t125 + t90) - g(2) * (-pkin(2) * t126 + t87) - g(3) * t128) * MDP(14) + t138 * t76 + t144 * (t129 * t143 - g(1) * t90 - g(2) * t87 - g(3) * (t107 * qJ(4) + t128)) + ((-pkin(4) * t101 + qJ(5) * t102) * MDP(22) - t113 * MDP(28) + t112 * MDP(29) - t140 * t101 - t139 * t102 + t141) * t77; (-MDP(14) - t144) * t76; -t144 * t77; (-g(1) * t78 + g(2) * t80 - t102 * t130) * MDP(22); (-g(1) * t71 - g(2) * t116) * MDP(28) + (g(1) * t72 - g(2) * t115) * MDP(29) + (-MDP(28) * t112 - MDP(29) * t113) * t130;];
taug  = t1;
