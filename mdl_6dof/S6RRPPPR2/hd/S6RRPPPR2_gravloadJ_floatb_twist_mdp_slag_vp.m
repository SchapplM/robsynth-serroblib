% Calculate Gravitation load on the joints for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPPPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:12:07
% EndTime: 2019-03-09 08:12:09
% DurationCPUTime: 0.53s
% Computational Cost: add. (239->82), mult. (289->109), div. (0->0), fcn. (250->10), ass. (0->49)
t84 = qJ(2) + pkin(9);
t78 = sin(t84);
t80 = cos(t84);
t118 = t80 * pkin(3) + t78 * qJ(4);
t117 = -MDP(14) + MDP(19);
t89 = sin(qJ(1));
t91 = cos(qJ(1));
t70 = g(1) * t91 + g(2) * t89;
t116 = t70 * t78;
t88 = sin(qJ(2));
t114 = pkin(2) * t88;
t111 = g(3) * t80;
t87 = -qJ(3) - pkin(7);
t110 = pkin(4) - t87;
t83 = pkin(10) + qJ(6);
t77 = sin(t83);
t109 = t77 * t91;
t79 = cos(t83);
t108 = t79 * t91;
t85 = sin(pkin(10));
t107 = t85 * t91;
t86 = cos(pkin(10));
t106 = t86 * t91;
t105 = t87 * t91;
t104 = t89 * t77;
t103 = t89 * t79;
t102 = t89 * t85;
t101 = t89 * t86;
t99 = qJ(4) * t80;
t98 = qJ(5) * t80;
t97 = -MDP(16) - MDP(20);
t90 = cos(qJ(2));
t81 = t90 * pkin(2);
t96 = t81 + t118;
t76 = t81 + pkin(1);
t72 = t91 * t76;
t95 = g(2) * (t118 * t91 + t72);
t94 = -pkin(3) * t78 - t114;
t69 = g(1) * t89 - g(2) * t91;
t93 = -t76 - t118;
t60 = -g(3) * t78 - t70 * t80;
t68 = t91 * t99;
t66 = t89 * t99;
t64 = -t78 * t104 + t108;
t63 = t78 * t103 + t109;
t62 = t78 * t109 + t103;
t61 = t78 * t108 - t104;
t59 = -t111 + t116;
t1 = [(-g(1) * (-t89 * t76 - t105) - g(2) * (-t89 * t87 + t72)) * MDP(12) + (g(1) * t105 - t95 + (-g(1) * t93 + g(2) * t87) * t89) * MDP(16) + (-g(1) * (-t78 * t102 + t106) - g(2) * (t78 * t107 + t101)) * MDP(17) + (-g(1) * (-t78 * t101 - t107) - g(2) * (t78 * t106 - t102)) * MDP(18) + (-t95 + (-g(1) * t110 - g(2) * t98) * t91 + (-g(1) * (t93 - t98) - g(2) * t110) * t89) * MDP(20) + (-g(1) * t64 - g(2) * t62) * MDP(26) + (g(1) * t63 - g(2) * t61) * MDP(27) + (MDP(3) - MDP(11) - MDP(13)) * t70 + (-t88 * MDP(10) + t78 * MDP(15) + t90 * MDP(9) + t117 * t80 + MDP(2)) * t69; (g(3) * t88 + t70 * t90) * MDP(10) + (-g(1) * (t94 * t91 + t68) - g(2) * (t94 * t89 + t66) - g(3) * t96) * MDP(16) + (-g(1) * (-t91 * t114 + t68) - g(2) * (-t89 * t114 + t66) - g(3) * (t96 + t98) + (pkin(3) + qJ(5)) * t116) * MDP(20) + (pkin(2) * MDP(12) + MDP(9)) * (-g(3) * t90 + t70 * t88) + t117 * t59 + (MDP(17) * t85 + MDP(18) * t86 + MDP(26) * t77 + MDP(27) * t79 + MDP(15)) * t60; (-MDP(12) + t97) * t69; t97 * t59; t60 * MDP(20); (-g(1) * t61 - g(2) * t63 + t79 * t111) * MDP(26) + (g(1) * t62 - g(2) * t64 - t77 * t111) * MDP(27);];
taug  = t1;
