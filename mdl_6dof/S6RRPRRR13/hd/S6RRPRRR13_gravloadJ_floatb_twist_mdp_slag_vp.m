% Calculate Gravitation load on the joints for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR13_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR13_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:53:20
% EndTime: 2019-03-09 14:53:22
% DurationCPUTime: 0.58s
% Computational Cost: add. (341->110), mult. (735->188), div. (0->0), fcn. (880->12), ass. (0->48)
t125 = MDP(9) - MDP(12);
t124 = -MDP(10) + MDP(13);
t103 = cos(pkin(6));
t94 = cos(qJ(1));
t101 = t94 * t103;
t89 = sin(qJ(2));
t90 = sin(qJ(1));
t93 = cos(qJ(2));
t78 = t89 * t101 + t90 * t93;
t87 = sin(qJ(5));
t91 = cos(qJ(5));
t86 = sin(pkin(6));
t109 = t86 * t94;
t77 = -t93 * t101 + t90 * t89;
t88 = sin(qJ(4));
t92 = cos(qJ(4));
t97 = t92 * t109 - t77 * t88;
t123 = t78 * t91 + t87 * t97;
t122 = -t78 * t87 + t91 * t97;
t85 = qJ(5) + qJ(6);
t83 = sin(t85);
t84 = cos(t85);
t121 = t78 * t84 + t83 * t97;
t120 = -t78 * t83 + t84 * t97;
t119 = g(3) * t86;
t114 = t83 * t88;
t113 = t84 * t88;
t112 = t86 * t89;
t111 = t86 * t90;
t110 = t86 * t93;
t108 = t87 * t88;
t107 = t88 * t89;
t106 = t88 * t91;
t105 = t89 * t91;
t102 = t90 * t103;
t79 = t93 * t102 + t94 * t89;
t70 = t92 * t111 + t79 * t88;
t80 = -t89 * t102 + t94 * t93;
t64 = -t70 * t83 + t80 * t84;
t65 = t70 * t84 + t80 * t83;
t76 = t103 * t92 - t88 * t110;
t104 = (-g(1) * t64 - g(2) * t121 - g(3) * (t84 * t112 - t76 * t83)) * MDP(34) + (g(1) * t65 - g(2) * t120 - g(3) * (-t83 * t112 - t76 * t84)) * MDP(35);
t71 = t88 * t109 + t77 * t92;
t68 = -g(1) * t79 - g(2) * t77 + g(3) * t110;
t69 = -t88 * t111 + t79 * t92;
t67 = t70 * t91 + t80 * t87;
t66 = -t70 * t87 + t80 * t91;
t1 = [(g(1) * t90 - g(2) * t94) * MDP(2) + (-g(1) * (-t90 * pkin(1) - t78 * pkin(2) + pkin(8) * t109 - t77 * qJ(3)) - g(2) * (pkin(1) * t94 + t80 * pkin(2) + pkin(8) * t111 + t79 * qJ(3))) * MDP(14) + (-g(1) * t97 - g(2) * t70) * MDP(20) + (g(1) * t71 - g(2) * t69) * MDP(21) + (-g(1) * t122 - g(2) * t67) * MDP(27) + (g(1) * t123 - g(2) * t66) * MDP(28) + (-g(1) * t120 - g(2) * t65) * MDP(34) + (g(1) * t121 - g(2) * t64) * MDP(35) + t125 * (g(1) * t78 - g(2) * t80) + (-t86 * MDP(11) + MDP(3)) * (g(1) * t94 + g(2) * t90) + t124 * (g(1) * t77 - g(2) * t79); (-g(1) * (-pkin(2) * t79 + qJ(3) * t80) - g(2) * (-pkin(2) * t77 + qJ(3) * t78) - (pkin(2) * t93 + qJ(3) * t89) * t119) * MDP(14) + (-g(1) * (t80 * t106 - t79 * t87) - g(2) * (t78 * t106 - t77 * t87) - (t88 * t105 + t87 * t93) * t119) * MDP(27) + (-g(1) * (-t80 * t108 - t79 * t91) - g(2) * (-t78 * t108 - t77 * t91) - (-t87 * t107 + t91 * t93) * t119) * MDP(28) + (-g(1) * (t80 * t113 - t79 * t83) - g(2) * (t78 * t113 - t77 * t83) - (t84 * t107 + t83 * t93) * t119) * MDP(34) + (-g(1) * (-t80 * t114 - t79 * t84) - g(2) * (-t78 * t114 - t77 * t84) - (-t83 * t107 + t84 * t93) * t119) * MDP(35) - t125 * t68 + (-t88 * MDP(20) - t92 * MDP(21) - t124) * (g(1) * t80 + g(2) * t78 + g(3) * t112); t68 * MDP(14); (g(1) * t70 - g(2) * t97 + g(3) * t76) * MDP(21) + (-MDP(27) * t91 + MDP(28) * t87 - MDP(34) * t84 + MDP(35) * t83 - MDP(20)) * (g(1) * t69 + g(2) * t71 + g(3) * (-t103 * t88 - t92 * t110)); (-g(1) * t66 - g(2) * t123 - g(3) * (t86 * t105 - t76 * t87)) * MDP(27) + (g(1) * t67 - g(2) * t122 - g(3) * (-t87 * t112 - t76 * t91)) * MDP(28) + t104; t104;];
taug  = t1;
