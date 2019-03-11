% Calculate Gravitation load on the joints for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:14:00
% EndTime: 2019-03-09 14:14:02
% DurationCPUTime: 0.61s
% Computational Cost: add. (466->107), mult. (670->180), div. (0->0), fcn. (784->14), ass. (0->47)
t91 = sin(pkin(6));
t98 = cos(qJ(1));
t112 = t91 * t98;
t109 = cos(pkin(6));
t106 = t98 * t109;
t94 = sin(qJ(2));
t95 = sin(qJ(1));
t97 = cos(qJ(2));
t77 = t94 * t106 + t95 * t97;
t89 = pkin(12) + qJ(4);
t88 = qJ(5) + t89;
t84 = sin(t88);
t85 = cos(t88);
t68 = -t84 * t112 + t77 * t85;
t76 = -t97 * t106 + t94 * t95;
t93 = sin(qJ(6));
t96 = cos(qJ(6));
t122 = t68 * t93 - t76 * t96;
t121 = t68 * t96 + t76 * t93;
t120 = MDP(10) - MDP(13);
t119 = g(3) * t91;
t116 = t85 * t93;
t115 = t85 * t96;
t114 = t91 * t94;
t113 = t91 * t95;
t111 = t93 * t97;
t110 = t96 * t97;
t86 = sin(t89);
t87 = cos(t89);
t108 = -t86 * t112 + t77 * t87;
t107 = t95 * t109;
t103 = t85 * t112 + t77 * t84;
t79 = -t94 * t107 + t97 * t98;
t70 = t85 * t113 - t79 * t84;
t71 = t84 * t113 + t79 * t85;
t75 = t109 * t84 + t85 * t114;
t105 = (g(1) * t71 + g(2) * t68 + g(3) * t75) * MDP(28) + (-MDP(34) * t96 + MDP(35) * t93 - MDP(27)) * (g(1) * t70 - g(2) * t103 + g(3) * (t109 * t85 - t84 * t114));
t102 = t87 * t112 + t77 * t86;
t78 = t97 * t107 + t98 * t94;
t100 = -g(1) * t78 - g(2) * t76 + t97 * t119;
t92 = cos(pkin(12));
t90 = sin(pkin(12));
t73 = t86 * t113 + t79 * t87;
t72 = t87 * t113 - t79 * t86;
t66 = t71 * t96 + t78 * t93;
t65 = -t71 * t93 + t78 * t96;
t1 = [(g(1) * t95 - g(2) * t98) * MDP(2) + (g(1) * t98 + g(2) * t95) * MDP(3) + (g(1) * t77 - g(2) * t79) * MDP(9) + (-g(1) * (t90 * t112 - t77 * t92) - g(2) * (t90 * t113 + t79 * t92)) * MDP(11) + (-g(1) * (t92 * t112 + t77 * t90) - g(2) * (t92 * t113 - t79 * t90)) * MDP(12) + (-g(1) * (-t95 * pkin(1) - t77 * pkin(2) + pkin(8) * t112 - t76 * qJ(3)) - g(2) * (pkin(1) * t98 + t79 * pkin(2) + pkin(8) * t113 + t78 * qJ(3))) * MDP(14) + (g(1) * t108 - g(2) * t73) * MDP(20) + (-g(1) * t102 - g(2) * t72) * MDP(21) + (g(1) * t68 - g(2) * t71) * MDP(27) + (-g(1) * t103 - g(2) * t70) * MDP(28) + (g(1) * t121 - g(2) * t66) * MDP(34) + (-g(1) * t122 - g(2) * t65) * MDP(35) - t120 * (g(1) * t76 - g(2) * t78); (-g(1) * (-pkin(2) * t78 + qJ(3) * t79) - g(2) * (-pkin(2) * t76 + qJ(3) * t77) - (pkin(2) * t97 + qJ(3) * t94) * t119) * MDP(14) + (-g(1) * (-t78 * t115 + t79 * t93) - g(2) * (-t76 * t115 + t77 * t93) - (t85 * t110 + t93 * t94) * t119) * MDP(34) + (-g(1) * (t78 * t116 + t79 * t96) - g(2) * (t76 * t116 + t77 * t96) - (-t85 * t111 + t94 * t96) * t119) * MDP(35) + t120 * (g(1) * t79 + g(2) * t77 + g(3) * t114) + (-t92 * MDP(11) + MDP(12) * t90 - MDP(20) * t87 + MDP(21) * t86 - t85 * MDP(27) + MDP(28) * t84 - MDP(9)) * t100; t100 * MDP(14); (-g(1) * t72 + g(2) * t102 - g(3) * (t109 * t87 - t86 * t114)) * MDP(20) + (g(1) * t73 + g(2) * t108 - g(3) * (-t109 * t86 - t87 * t114)) * MDP(21) + t105; t105; (-g(1) * t65 + g(2) * t122 - g(3) * (-t91 * t110 - t75 * t93)) * MDP(34) + (g(1) * t66 + g(2) * t121 - g(3) * (t91 * t111 - t75 * t96)) * MDP(35);];
taug  = t1;
