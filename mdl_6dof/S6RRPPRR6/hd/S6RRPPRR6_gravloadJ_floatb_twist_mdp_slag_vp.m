% Calculate Gravitation load on the joints for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:15:16
% EndTime: 2019-03-09 09:15:17
% DurationCPUTime: 0.47s
% Computational Cost: add. (245->81), mult. (404->123), div. (0->0), fcn. (413->10), ass. (0->47)
t104 = sin(qJ(6));
t107 = cos(qJ(6));
t105 = sin(qJ(2));
t108 = cos(qJ(2));
t121 = pkin(10) + qJ(5);
t120 = cos(t121);
t95 = sin(t121);
t110 = t105 * t95 + t108 * t120;
t118 = t105 * t120;
t132 = g(3) * (-t108 * t95 + t118);
t106 = sin(qJ(1));
t125 = t106 * t108;
t71 = -t106 * t118 + t95 * t125;
t72 = t110 * t106;
t109 = cos(qJ(1));
t124 = t108 * t109;
t73 = -t109 * t118 + t95 * t124;
t74 = t110 * t109;
t139 = (MDP(31) * t107 - MDP(32) * t104 + MDP(24)) * (g(1) * t73 + g(2) * t71 + g(3) * t110) + (g(1) * t74 + g(2) * t72 + t132) * MDP(25);
t138 = MDP(9) + MDP(11);
t137 = MDP(10) - MDP(13);
t127 = t108 * pkin(2) + t105 * qJ(3);
t88 = g(1) * t109 + g(2) * t106;
t135 = t88 * t105;
t131 = pkin(3) * t108;
t130 = g(1) * t106;
t126 = t105 * t109;
t119 = t109 * pkin(1) + pkin(2) * t124 + t106 * pkin(7) + qJ(3) * t126;
t87 = -g(2) * t109 + t130;
t117 = t72 * t104 - t107 * t109;
t116 = t104 * t109 + t72 * t107;
t102 = sin(pkin(10));
t103 = cos(pkin(10));
t115 = t102 * t108 - t103 * t105;
t114 = t102 * t105 + t103 * t108;
t113 = -pkin(1) - t127;
t99 = t109 * pkin(7);
t93 = qJ(3) * t124;
t91 = qJ(3) * t125;
t80 = t114 * t109;
t79 = t115 * t109;
t78 = t114 * t106;
t77 = t115 * t106;
t75 = -g(3) * t108 + t135;
t70 = -t104 * t106 + t107 * t74;
t69 = -t104 * t74 - t106 * t107;
t1 = [(-g(1) * t99 - g(2) * t119 - t113 * t130) * MDP(14) + (g(1) * t78 - g(2) * t80) * MDP(15) + (-g(1) * t77 + g(2) * t79) * MDP(16) + (-g(1) * (-qJ(4) * t109 + t99) - g(2) * (pkin(3) * t124 + t119) + (-g(1) * (t113 - t131) + g(2) * qJ(4)) * t106) * MDP(18) + (g(1) * t72 - g(2) * t74) * MDP(24) + (-g(1) * t71 + g(2) * t73) * MDP(25) + (g(1) * t116 - g(2) * t70) * MDP(31) + (-g(1) * t117 - g(2) * t69) * MDP(32) + (MDP(3) - MDP(12) + MDP(17)) * t88 + (-t137 * t105 + t138 * t108 + MDP(2)) * t87; (-g(1) * (-pkin(2) * t126 + t93) - g(2) * (-pkin(2) * t105 * t106 + t91) - g(3) * t127) * MDP(14) + (-g(1) * t79 - g(2) * t77 - g(3) * t114) * MDP(15) + (-g(1) * t80 - g(2) * t78 + g(3) * t115) * MDP(16) + (-g(1) * t93 - g(2) * t91 - g(3) * (t127 + t131) + (pkin(2) + pkin(3)) * t135) * MDP(18) + t137 * (g(3) * t105 + t88 * t108) + t138 * t75 - t139; (-MDP(14) - MDP(18)) * t75; t87 * MDP(18); t139; (-g(1) * t69 + g(2) * t117 + t104 * t132) * MDP(31) + (g(1) * t70 + g(2) * t116 + t107 * t132) * MDP(32);];
taug  = t1;
