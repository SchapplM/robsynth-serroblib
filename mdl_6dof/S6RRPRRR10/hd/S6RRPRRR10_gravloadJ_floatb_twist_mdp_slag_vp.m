% Calculate Gravitation load on the joints for
% S6RRPRRR10
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
%   see S6RRPRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:25:35
% EndTime: 2019-03-09 14:25:37
% DurationCPUTime: 0.69s
% Computational Cost: add. (466->119), mult. (768->204), div. (0->0), fcn. (920->14), ass. (0->51)
t122 = cos(qJ(1));
t91 = sin(pkin(6));
t105 = t91 * t122;
t106 = cos(pkin(6));
t102 = t106 * t122;
t94 = sin(qJ(2));
t95 = sin(qJ(1));
t97 = cos(qJ(2));
t78 = t102 * t94 + t95 * t97;
t88 = pkin(12) + qJ(4);
t84 = sin(t88);
t85 = cos(t88);
t71 = -t84 * t105 + t78 * t85;
t77 = -t102 * t97 + t95 * t94;
t93 = sin(qJ(5));
t96 = cos(qJ(5));
t128 = t71 * t93 - t77 * t96;
t127 = t71 * t96 + t77 * t93;
t89 = qJ(5) + qJ(6);
t86 = sin(t89);
t87 = cos(t89);
t126 = t71 * t86 - t77 * t87;
t125 = t71 * t87 + t77 * t86;
t124 = MDP(10) - MDP(13);
t123 = g(3) * t91;
t117 = t85 * t86;
t116 = t85 * t87;
t115 = t85 * t93;
t114 = t85 * t96;
t113 = t85 * t97;
t112 = t91 * t94;
t111 = t91 * t95;
t110 = t91 * t97;
t109 = t93 * t97;
t108 = t96 * t97;
t104 = t95 * t106;
t80 = -t104 * t94 + t122 * t97;
t74 = t111 * t84 + t80 * t85;
t79 = t104 * t97 + t122 * t94;
t66 = -t74 * t86 + t79 * t87;
t67 = t74 * t87 + t79 * t86;
t76 = t106 * t84 + t112 * t85;
t107 = (-g(1) * t66 + g(2) * t126 - g(3) * (-t110 * t87 - t76 * t86)) * MDP(34) + (g(1) * t67 + g(2) * t125 - g(3) * (t110 * t86 - t76 * t87)) * MDP(35);
t100 = t105 * t85 + t78 * t84;
t99 = -g(1) * t79 - g(2) * t77 + g(3) * t110;
t92 = cos(pkin(12));
t90 = sin(pkin(12));
t73 = t111 * t85 - t80 * t84;
t69 = t74 * t96 + t79 * t93;
t68 = -t74 * t93 + t79 * t96;
t1 = [(g(1) * t95 - g(2) * t122) * MDP(2) + (g(1) * t122 + g(2) * t95) * MDP(3) + (g(1) * t78 - g(2) * t80) * MDP(9) + (-g(1) * (t105 * t90 - t78 * t92) - g(2) * (t111 * t90 + t80 * t92)) * MDP(11) + (-g(1) * (t105 * t92 + t78 * t90) - g(2) * (t111 * t92 - t80 * t90)) * MDP(12) + (-g(1) * (-t95 * pkin(1) - t78 * pkin(2) + pkin(8) * t105 - t77 * qJ(3)) - g(2) * (pkin(1) * t122 + t80 * pkin(2) + pkin(8) * t111 + t79 * qJ(3))) * MDP(14) + (g(1) * t71 - g(2) * t74) * MDP(20) + (-g(1) * t100 - g(2) * t73) * MDP(21) + (g(1) * t127 - g(2) * t69) * MDP(27) + (-g(1) * t128 - g(2) * t68) * MDP(28) + (g(1) * t125 - g(2) * t67) * MDP(34) + (-g(1) * t126 - g(2) * t66) * MDP(35) - t124 * (g(1) * t77 - g(2) * t79); (-g(1) * (-t79 * pkin(2) + t80 * qJ(3)) - g(2) * (-t77 * pkin(2) + t78 * qJ(3)) - (pkin(2) * t97 + qJ(3) * t94) * t123) * MDP(14) + (-g(1) * (-t114 * t79 + t80 * t93) - g(2) * (-t114 * t77 + t78 * t93) - (t108 * t85 + t93 * t94) * t123) * MDP(27) + (-g(1) * (t115 * t79 + t80 * t96) - g(2) * (t115 * t77 + t78 * t96) - (-t109 * t85 + t94 * t96) * t123) * MDP(28) + (-g(1) * (-t116 * t79 + t80 * t86) - g(2) * (-t116 * t77 + t78 * t86) - (t113 * t87 + t86 * t94) * t123) * MDP(34) + (-g(1) * (t117 * t79 + t80 * t87) - g(2) * (t117 * t77 + t78 * t87) - (-t113 * t86 + t87 * t94) * t123) * MDP(35) + t124 * (g(1) * t80 + g(2) * t78 + g(3) * t112) + (-t92 * MDP(11) + t90 * MDP(12) - MDP(20) * t85 + t84 * MDP(21) - MDP(9)) * t99; t99 * MDP(14); (g(1) * t74 + g(2) * t71 + g(3) * t76) * MDP(21) + (-MDP(27) * t96 + MDP(28) * t93 - MDP(34) * t87 + MDP(35) * t86 - MDP(20)) * (g(1) * t73 - g(2) * t100 + g(3) * (t106 * t85 - t112 * t84)); (-g(1) * t68 + g(2) * t128 - g(3) * (-t108 * t91 - t76 * t93)) * MDP(27) + (g(1) * t69 + g(2) * t127 - g(3) * (t109 * t91 - t76 * t96)) * MDP(28) + t107; t107;];
taug  = t1;
