% Calculate Gravitation load on the joints for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:21:00
% EndTime: 2019-03-09 09:21:01
% DurationCPUTime: 0.56s
% Computational Cost: add. (244->94), mult. (592->153), div. (0->0), fcn. (675->10), ass. (0->42)
t131 = -MDP(10) + MDP(13) + MDP(15);
t130 = MDP(16) - MDP(9) - MDP(11);
t100 = cos(qJ(5));
t102 = cos(qJ(1));
t94 = sin(pkin(6));
t116 = t102 * t94;
t101 = cos(qJ(2));
t115 = cos(pkin(6));
t111 = t102 * t115;
t97 = sin(qJ(2));
t98 = sin(qJ(1));
t82 = -t101 * t111 + t97 * t98;
t96 = sin(qJ(5));
t72 = t82 * t100 + t96 * t116;
t83 = t98 * t101 + t97 * t111;
t95 = sin(qJ(6));
t99 = cos(qJ(6));
t129 = t72 * t95 - t83 * t99;
t128 = t72 * t99 + t83 * t95;
t125 = g(3) * t94;
t122 = t94 * t97;
t121 = t94 * t98;
t117 = t101 * t94;
t120 = pkin(2) * t117 + qJ(3) * t122;
t119 = t100 * t95;
t118 = t100 * t99;
t114 = t98 * t115;
t113 = -t82 * pkin(2) + qJ(3) * t83;
t84 = t101 * t114 + t102 * t97;
t85 = t101 * t102 - t97 * t114;
t112 = -t84 * pkin(2) + qJ(3) * t85;
t108 = g(1) * t98 - g(2) * t102;
t106 = t102 * pkin(1) + t85 * pkin(2) + pkin(8) * t121 + qJ(3) * t84;
t71 = t100 * t116 - t82 * t96;
t105 = -t98 * pkin(1) - t83 * pkin(2) + pkin(8) * t116 - t82 * qJ(3);
t66 = -g(1) * t84 - g(2) * t82 + g(3) * t117;
t81 = t100 * t117 + t115 * t96;
t75 = t100 * t84 - t96 * t121;
t74 = -t100 * t121 - t84 * t96;
t65 = t75 * t99 + t85 * t95;
t64 = -t75 * t95 + t85 * t99;
t1 = [t108 * MDP(2) + (-g(1) * t105 - g(2) * t106) * MDP(14) + (-g(1) * (-t83 * pkin(3) - qJ(4) * t116 + t105) - g(2) * (pkin(3) * t85 - qJ(4) * t121 + t106)) * MDP(18) + (g(1) * t72 - g(2) * t75) * MDP(24) + (g(1) * t71 - g(2) * t74) * MDP(25) + (g(1) * t128 - g(2) * t65) * MDP(31) + (-g(1) * t129 - g(2) * t64) * MDP(32) + t131 * (g(1) * t82 - g(2) * t84) - t130 * (g(1) * t83 - g(2) * t85) + (MDP(3) + (-MDP(12) + MDP(17)) * t94) * (g(1) * t102 + g(2) * t98); (-g(1) * t112 - g(2) * t113 - g(3) * t120) * MDP(14) + (-g(1) * (-pkin(3) * t84 + t112) - g(2) * (-pkin(3) * t82 + t113) - g(3) * (pkin(3) * t117 + t120)) * MDP(18) + (-g(1) * (t85 * t118 - t84 * t95) - g(2) * (t83 * t118 - t82 * t95) - (t101 * t95 + t97 * t118) * t125) * MDP(31) + (-g(1) * (-t85 * t119 - t84 * t99) - g(2) * (-t83 * t119 - t82 * t99) - (t101 * t99 - t97 * t119) * t125) * MDP(32) + t130 * t66 + (-t100 * MDP(24) + MDP(25) * t96 - t131) * (g(1) * t85 + g(2) * t83 + g(3) * t122); (MDP(14) + MDP(18)) * t66; (g(3) * t115 + t108 * t94) * MDP(18); (g(1) * t75 + g(2) * t72 - g(3) * t81) * MDP(25) + (-MDP(31) * t99 + MDP(32) * t95 - MDP(24)) * (g(1) * t74 + g(2) * t71 + g(3) * (-t115 * t100 + t96 * t117)); (-g(1) * t64 + g(2) * t129 - g(3) * (t99 * t122 + t81 * t95)) * MDP(31) + (g(1) * t65 + g(2) * t128 - g(3) * (-t95 * t122 + t81 * t99)) * MDP(32);];
taug  = t1;
