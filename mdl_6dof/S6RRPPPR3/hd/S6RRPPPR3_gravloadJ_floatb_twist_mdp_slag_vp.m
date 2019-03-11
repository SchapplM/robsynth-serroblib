% Calculate Gravitation load on the joints for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:15:39
% EndTime: 2019-03-09 08:15:41
% DurationCPUTime: 0.62s
% Computational Cost: add. (181->85), mult. (340->114), div. (0->0), fcn. (297->8), ass. (0->50)
t90 = sin(qJ(2));
t80 = t90 * qJ(3);
t123 = pkin(1) + t80;
t122 = -MDP(10) + MDP(13) + MDP(15);
t121 = MDP(9) + MDP(11) - MDP(16) + MDP(21);
t91 = sin(qJ(1));
t116 = g(2) * t91;
t93 = cos(qJ(1));
t72 = g(1) * t93 + t116;
t120 = t72 * t90;
t92 = cos(qJ(2));
t68 = g(3) * t90 + t72 * t92;
t119 = pkin(2) + pkin(3);
t118 = g(1) * t91;
t114 = g(3) * t92;
t83 = t92 * pkin(2);
t82 = t92 * pkin(3);
t113 = g(2) * qJ(4);
t112 = t90 * t93;
t87 = pkin(9) + qJ(6);
t78 = sin(t87);
t111 = t91 * t78;
t79 = cos(t87);
t110 = t91 * t79;
t88 = sin(pkin(9));
t109 = t91 * t88;
t89 = cos(pkin(9));
t108 = t91 * t89;
t107 = t92 * t93;
t106 = t83 + t80;
t105 = qJ(3) * t92;
t104 = qJ(5) * t92;
t103 = MDP(18) + MDP(22);
t102 = qJ(5) + t119;
t101 = t82 + t106;
t99 = pkin(2) * t107 + t91 * pkin(7) + t123 * t93;
t98 = g(1) * t102;
t97 = pkin(3) * t107 + t99;
t71 = -g(2) * t93 + t118;
t84 = t93 * pkin(7);
t96 = g(1) * (-qJ(4) * t93 + t84);
t95 = -t123 - t83;
t75 = t93 * t105;
t73 = t91 * t105;
t67 = -t114 + t120;
t66 = t112 * t79 - t111;
t65 = -t112 * t78 - t110;
t64 = -t110 * t90 - t78 * t93;
t63 = t111 * t90 - t79 * t93;
t1 = [(-g(1) * t84 - g(2) * t99 - t118 * t95) * MDP(14) + (-t96 - g(2) * t97 + (-g(1) * (t95 - t82) + t113) * t91) * MDP(18) + (-g(1) * (-t108 * t90 - t88 * t93) - g(2) * (t112 * t89 - t109)) * MDP(19) + (-g(1) * (t109 * t90 - t89 * t93) - g(2) * (-t112 * t88 - t108)) * MDP(20) + (-t96 - g(2) * (pkin(4) * t112 + t104 * t93 + t97) + (-g(1) * (-pkin(4) * t90 - t123) + t113 + t92 * t98) * t91) * MDP(22) + (-g(1) * t64 - g(2) * t66) * MDP(28) + (-g(1) * t63 - g(2) * t65) * MDP(29) + (MDP(3) - MDP(12) + MDP(17)) * t72 + (t121 * t92 + t122 * t90 + MDP(2)) * t71; (-g(1) * (-pkin(2) * t112 + t75) - g(2) * (-pkin(2) * t90 * t91 + t73) - g(3) * t106) * MDP(14) + (-g(1) * t75 - g(2) * t73 - g(3) * t101 + t119 * t120) * MDP(18) + (-g(1) * (pkin(4) * t107 + t75) - g(2) * (pkin(4) * t91 * t92 + t73) - g(3) * (t101 + t104) + (-g(3) * pkin(4) + t102 * t116 + t93 * t98) * t90) * MDP(22) + t121 * t67 + (-MDP(19) * t89 + MDP(20) * t88 - MDP(28) * t79 + MDP(29) * t78 - t122) * t68; (-MDP(14) - t103) * t67; t103 * t71; -t68 * MDP(22); (-g(1) * t65 + g(2) * t63 - t114 * t78) * MDP(28) + (g(1) * t66 - g(2) * t64 - t114 * t79) * MDP(29);];
taug  = t1;
