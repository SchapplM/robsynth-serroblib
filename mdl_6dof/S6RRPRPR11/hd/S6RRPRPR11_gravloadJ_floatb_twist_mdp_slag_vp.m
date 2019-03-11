% Calculate Gravitation load on the joints for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:15:04
% EndTime: 2019-03-09 11:15:05
% DurationCPUTime: 0.42s
% Computational Cost: add. (212->77), mult. (321->109), div. (0->0), fcn. (293->8), ass. (0->45)
t130 = MDP(10) - MDP(13);
t129 = MDP(9) - MDP(12) + MDP(22);
t96 = sin(qJ(2));
t99 = cos(qJ(2));
t113 = t99 * pkin(2) + t96 * qJ(3);
t100 = cos(qJ(1));
t97 = sin(qJ(1));
t128 = -g(1) * t100 - g(2) * t97;
t72 = g(3) * t96 - t128 * t99;
t125 = g(1) * t97;
t121 = g(3) * t99;
t95 = sin(qJ(4));
t119 = t95 * t99;
t87 = qJ(4) + pkin(10) + qJ(6);
t84 = sin(t87);
t118 = t97 * t84;
t85 = cos(t87);
t117 = t97 * t85;
t116 = t97 * t95;
t98 = cos(qJ(4));
t115 = t97 * t98;
t112 = t100 * t96;
t67 = t85 * t112 - t118;
t68 = t84 * t112 + t117;
t69 = t100 * t84 + t96 * t117;
t70 = t100 * t85 - t96 * t118;
t114 = (-g(1) * t67 - g(2) * t69 + t85 * t121) * MDP(29) + (g(1) * t68 - g(2) * t70 - t84 * t121) * MDP(30);
t111 = t100 * t98;
t110 = t100 * t99;
t109 = qJ(3) * t100;
t107 = t95 * t112;
t106 = t100 * pkin(1) + pkin(2) * t110 + t97 * pkin(7) + t96 * t109;
t94 = -qJ(5) - pkin(8);
t104 = pkin(4) * t95 * t96 - t94 * t99;
t103 = -pkin(1) - t113;
t75 = t100 * t95 + t96 * t115;
t73 = t96 * t111 - t116;
t91 = t100 * pkin(7);
t86 = pkin(4) * t98 + pkin(3);
t82 = t99 * t109;
t80 = t97 * t99 * qJ(3);
t76 = -t96 * t116 + t111;
t74 = t107 + t115;
t71 = -t128 * t96 - t121;
t1 = [(-g(1) * t91 - g(2) * t106 - t103 * t125) * MDP(14) + (-g(1) * t76 - g(2) * t74) * MDP(20) + (g(1) * t75 - g(2) * t73) * MDP(21) + (-g(1) * (t100 * t86 + t91) - g(2) * (pkin(4) * t107 - t94 * t110 + t106) + (-g(1) * (t103 - t104) - g(2) * t86) * t97) * MDP(23) + (-g(1) * t70 - g(2) * t68) * MDP(29) + (g(1) * t69 - g(2) * t67) * MDP(30) - (MDP(3) - MDP(11)) * t128 + (t129 * t99 - t130 * t96 + MDP(2)) * (-g(2) * t100 + t125); (-g(1) * (-pkin(2) * t112 + t82) - g(2) * (-pkin(2) * t96 * t97 + t80) - g(3) * t113) * MDP(14) + (-g(1) * t82 - g(2) * t80 - g(3) * (t104 + t113) + t128 * (pkin(4) * t119 + (-pkin(2) + t94) * t96)) * MDP(23) + t129 * t71 + (-t95 * MDP(20) - MDP(21) * t98 - MDP(29) * t84 - MDP(30) * t85 + t130) * t72; (-MDP(14) - MDP(23)) * t71; (g(1) * t74 - g(2) * t76 - g(3) * t119) * MDP(21) + t114 + (pkin(4) * MDP(23) + MDP(20)) * (-g(1) * t73 - g(2) * t75 + t98 * t121); -t72 * MDP(23); t114;];
taug  = t1;
