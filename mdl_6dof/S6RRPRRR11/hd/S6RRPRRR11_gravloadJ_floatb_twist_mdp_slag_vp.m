% Calculate Gravitation load on the joints for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:32:50
% EndTime: 2019-03-09 14:32:51
% DurationCPUTime: 0.27s
% Computational Cost: add. (261->70), mult. (328->101), div. (0->0), fcn. (321->10), ass. (0->46)
t86 = sin(qJ(1));
t89 = cos(qJ(1));
t76 = g(1) * t89 + g(2) * t86;
t115 = MDP(9) - MDP(12);
t114 = MDP(10) - MDP(13);
t88 = cos(qJ(2));
t109 = g(3) * t88;
t83 = qJ(4) + qJ(5);
t81 = qJ(6) + t83;
t77 = sin(t81);
t108 = t86 * t77;
t78 = cos(t81);
t107 = t86 * t78;
t79 = sin(t83);
t106 = t86 * t79;
t80 = cos(t83);
t105 = t86 * t80;
t84 = sin(qJ(4));
t104 = t86 * t84;
t87 = cos(qJ(4));
t103 = t86 * t87;
t102 = t89 * t77;
t101 = t89 * t78;
t100 = t89 * t79;
t99 = t89 * t80;
t98 = t89 * t84;
t97 = t89 * t87;
t85 = sin(qJ(2));
t60 = t85 * t101 - t108;
t61 = t85 * t102 + t107;
t62 = t85 * t107 + t102;
t63 = -t85 * t108 + t101;
t96 = (-g(1) * t60 - g(2) * t62 + t78 * t109) * MDP(34) + (g(1) * t61 - g(2) * t63 - t77 * t109) * MDP(35);
t64 = t85 * t99 - t106;
t65 = t85 * t100 + t105;
t66 = t85 * t105 + t100;
t67 = -t85 * t106 + t99;
t95 = (-g(1) * t64 - g(2) * t66 + t80 * t109) * MDP(27) + (g(1) * t65 - g(2) * t67 - t79 * t109) * MDP(28) + t96;
t93 = t88 * pkin(2) + t85 * qJ(3);
t91 = pkin(1) + t93;
t73 = -t85 * t104 + t97;
t72 = t85 * t103 + t98;
t71 = t85 * t98 + t103;
t70 = t85 * t97 - t104;
t68 = t76 * t85 - t109;
t1 = [((-g(1) * pkin(7) - g(2) * t91) * t89 + (-g(2) * pkin(7) + g(1) * t91) * t86) * MDP(14) + (-g(1) * t73 - g(2) * t71) * MDP(20) + (g(1) * t72 - g(2) * t70) * MDP(21) + (-g(1) * t67 - g(2) * t65) * MDP(27) + (g(1) * t66 - g(2) * t64) * MDP(28) + (-g(1) * t63 - g(2) * t61) * MDP(34) + (g(1) * t62 - g(2) * t60) * MDP(35) + (MDP(3) - MDP(11)) * t76 + (-t114 * t85 + t115 * t88 + MDP(2)) * (g(1) * t86 - g(2) * t89); (-g(3) * t93 + t76 * (pkin(2) * t85 - qJ(3) * t88)) * MDP(14) + t115 * t68 + (-t84 * MDP(20) - t87 * MDP(21) - t79 * MDP(27) - t80 * MDP(28) - t77 * MDP(34) - t78 * MDP(35) + t114) * (g(3) * t85 + t76 * t88); -t68 * MDP(14); (-g(1) * t70 - g(2) * t72 + t87 * t109) * MDP(20) + (g(1) * t71 - g(2) * t73 - t84 * t109) * MDP(21) + t95; t95; t96;];
taug  = t1;
