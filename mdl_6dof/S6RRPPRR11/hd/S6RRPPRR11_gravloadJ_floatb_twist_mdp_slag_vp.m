% Calculate Gravitation load on the joints for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:42:45
% EndTime: 2019-03-09 09:42:47
% DurationCPUTime: 0.73s
% Computational Cost: add. (321->104), mult. (630->169), div. (0->0), fcn. (722->12), ass. (0->47)
t138 = -MDP(10) + MDP(13);
t137 = MDP(12) - MDP(9) - MDP(17);
t103 = sin(qJ(6));
t106 = cos(qJ(6));
t101 = sin(pkin(6));
t108 = cos(qJ(1));
t123 = t101 * t108;
t104 = sin(qJ(2));
t105 = sin(qJ(1));
t107 = cos(qJ(2));
t127 = cos(pkin(6));
t117 = t108 * t127;
t86 = t104 * t105 - t107 * t117;
t99 = pkin(11) + qJ(5);
t96 = sin(t99);
t97 = cos(t99);
t112 = t97 * t123 - t86 * t96;
t87 = t104 * t117 + t105 * t107;
t135 = t103 * t112 + t106 * t87;
t134 = -t103 * t87 + t106 * t112;
t133 = g(3) * t101;
t124 = t101 * t107;
t126 = t101 * t104;
t132 = pkin(2) * t124 + qJ(3) * t126;
t130 = t103 * t96;
t128 = t106 * t96;
t125 = t101 * t105;
t122 = t103 * t104;
t121 = t104 * t106;
t120 = -t86 * pkin(2) + qJ(3) * t87;
t118 = t105 * t127;
t88 = t108 * t104 + t107 * t118;
t89 = -t104 * t118 + t107 * t108;
t119 = -t88 * pkin(2) + qJ(3) * t89;
t113 = t108 * pkin(1) + t89 * pkin(2) + pkin(8) * t125 + qJ(3) * t88;
t76 = t96 * t123 + t86 * t97;
t110 = -t105 * pkin(1) - t87 * pkin(2) + pkin(8) * t123 - t86 * qJ(3);
t109 = g(1) * t89 + g(2) * t87 + g(3) * t126;
t70 = -g(1) * t88 - g(2) * t86 + g(3) * t124;
t102 = cos(pkin(11));
t100 = sin(pkin(11));
t81 = -t96 * t124 + t127 * t97;
t75 = t97 * t125 + t88 * t96;
t74 = -t96 * t125 + t88 * t97;
t69 = t103 * t89 + t106 * t75;
t68 = -t103 * t75 + t106 * t89;
t1 = [(g(1) * t105 - g(2) * t108) * MDP(2) + (-g(1) * t110 - g(2) * t113) * MDP(14) + (-g(1) * (-t86 * t100 + t102 * t123) - g(2) * (t100 * t88 + t102 * t125)) * MDP(15) + (-g(1) * (-t100 * t123 - t86 * t102) - g(2) * (-t100 * t125 + t102 * t88)) * MDP(16) + (-g(1) * (pkin(3) * t123 - t87 * qJ(4) + t110) - g(2) * (pkin(3) * t125 + qJ(4) * t89 + t113)) * MDP(18) + (-g(1) * t112 - g(2) * t75) * MDP(24) + (g(1) * t76 - g(2) * t74) * MDP(25) + (-g(1) * t134 - g(2) * t69) * MDP(31) + (g(1) * t135 - g(2) * t68) * MDP(32) + t138 * (g(1) * t86 - g(2) * t88) - t137 * (g(1) * t87 - g(2) * t89) + (-t101 * MDP(11) + MDP(3)) * (g(1) * t108 + g(2) * t105); (-g(1) * t119 - g(2) * t120 - g(3) * t132) * MDP(14) + (-g(1) * (-qJ(4) * t88 + t119) - g(2) * (-qJ(4) * t86 + t120) - g(3) * (qJ(4) * t124 + t132)) * MDP(18) + (-g(1) * (-t103 * t88 + t89 * t128) - g(2) * (-t103 * t86 + t87 * t128) - (t103 * t107 + t96 * t121) * t133) * MDP(31) + (-g(1) * (-t106 * t88 - t89 * t130) - g(2) * (-t106 * t86 - t87 * t130) - (t106 * t107 - t96 * t122) * t133) * MDP(32) + t137 * t70 + (-MDP(15) * t100 - MDP(16) * t102 - t96 * MDP(24) - MDP(25) * t97 - t138) * t109; (MDP(14) + MDP(18)) * t70; -t109 * MDP(18); (g(1) * t75 - g(2) * t112 + g(3) * t81) * MDP(25) + (-MDP(31) * t106 + MDP(32) * t103 - MDP(24)) * (g(1) * t74 + g(2) * t76 + g(3) * (-t97 * t124 - t127 * t96)); (-g(1) * t68 - g(2) * t135 - g(3) * (t101 * t121 - t103 * t81)) * MDP(31) + (g(1) * t69 - g(2) * t134 - g(3) * (-t101 * t122 - t106 * t81)) * MDP(32);];
taug  = t1;
