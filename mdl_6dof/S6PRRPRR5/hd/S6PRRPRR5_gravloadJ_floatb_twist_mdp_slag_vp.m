% Calculate Gravitation load on the joints for
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:20:48
% EndTime: 2019-03-08 22:20:50
% DurationCPUTime: 0.73s
% Computational Cost: add. (389->104), mult. (718->177), div. (0->0), fcn. (861->14), ass. (0->46)
t114 = MDP(11) - MDP(14);
t85 = sin(qJ(2));
t87 = cos(qJ(2));
t98 = cos(pkin(11));
t99 = cos(pkin(6));
t93 = t99 * t98;
t97 = sin(pkin(11));
t64 = t85 * t97 - t87 * t93;
t92 = t99 * t97;
t66 = t85 * t98 + t87 * t92;
t113 = -g(1) * t66 - g(2) * t64;
t82 = sin(pkin(6));
t110 = g(3) * t82;
t80 = pkin(12) + qJ(5);
t79 = qJ(6) + t80;
t75 = sin(t79);
t86 = cos(qJ(3));
t109 = t75 * t86;
t76 = cos(t79);
t108 = t76 * t86;
t77 = sin(t80);
t107 = t77 * t86;
t78 = cos(t80);
t106 = t78 * t86;
t81 = sin(pkin(12));
t105 = t81 * t86;
t104 = t82 * t85;
t103 = t82 * t87;
t83 = cos(pkin(12));
t102 = t83 * t86;
t101 = t86 * t87;
t65 = t85 * t93 + t87 * t97;
t84 = sin(qJ(3));
t96 = t82 * t98;
t61 = t65 * t86 - t84 * t96;
t67 = -t85 * t92 + t87 * t98;
t95 = t82 * t97;
t63 = t67 * t86 + t84 * t95;
t69 = t104 * t86 + t84 * t99;
t100 = (-g(1) * (-t63 * t75 + t66 * t76) - g(2) * (-t61 * t75 + t64 * t76) - g(3) * (-t103 * t76 - t69 * t75)) * MDP(28) + (-g(1) * (-t63 * t76 - t66 * t75) - g(2) * (-t61 * t76 - t64 * t75) - g(3) * (t103 * t75 - t69 * t76)) * MDP(29);
t94 = -g(1) * t67 - g(2) * t65;
t60 = t65 * t84 + t86 * t96;
t62 = t67 * t84 - t86 * t95;
t68 = t104 * t84 - t86 * t99;
t90 = g(1) * t62 + g(2) * t60 + g(3) * t68;
t1 = [(-MDP(1) - MDP(15)) * g(3); (g(3) * t104 - t94) * MDP(4) + (-g(1) * (-t102 * t66 + t67 * t81) - g(2) * (-t102 * t64 + t65 * t81) - (t101 * t83 + t81 * t85) * t110) * MDP(12) + (-g(1) * (t105 * t66 + t67 * t83) - g(2) * (t105 * t64 + t65 * t83) - (-t101 * t81 + t83 * t85) * t110) * MDP(13) + ((-t110 * t85 + t94) * pkin(8) + (-t110 * t87 - t113) * (pkin(3) * t86 + qJ(4) * t84 + pkin(2))) * MDP(15) + (-g(1) * (-t106 * t66 + t67 * t77) - g(2) * (-t106 * t64 + t65 * t77) - (t101 * t78 + t77 * t85) * t110) * MDP(21) + (-g(1) * (t107 * t66 + t67 * t78) - g(2) * (t107 * t64 + t65 * t78) - (-t101 * t77 + t78 * t85) * t110) * MDP(22) + (-g(1) * (-t108 * t66 + t67 * t75) - g(2) * (-t108 * t64 + t65 * t75) - (t101 * t76 + t75 * t85) * t110) * MDP(28) + (-g(1) * (t109 * t66 + t67 * t76) - g(2) * (t109 * t64 + t65 * t76) - (-t101 * t75 + t76 * t85) * t110) * MDP(29) + (-t86 * MDP(10) + t114 * t84 - MDP(3)) * (g(3) * t103 + t113); (-g(1) * (-pkin(3) * t62 + qJ(4) * t63) - g(2) * (-pkin(3) * t60 + qJ(4) * t61) - g(3) * (-pkin(3) * t68 + qJ(4) * t69)) * MDP(15) + t114 * (g(1) * t63 + g(2) * t61 + g(3) * t69) + (MDP(12) * t83 - MDP(13) * t81 + MDP(21) * t78 - MDP(22) * t77 + MDP(28) * t76 - MDP(29) * t75 + MDP(10)) * t90; -t90 * MDP(15); (-g(1) * (-t63 * t77 + t66 * t78) - g(2) * (-t61 * t77 + t64 * t78) - g(3) * (-t103 * t78 - t69 * t77)) * MDP(21) + (-g(1) * (-t63 * t78 - t66 * t77) - g(2) * (-t61 * t78 - t64 * t77) - g(3) * (t103 * t77 - t69 * t78)) * MDP(22) + t100; t100;];
taug  = t1;
