% Calculate Gravitation load on the joints for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:34:12
% EndTime: 2019-03-08 22:34:14
% DurationCPUTime: 0.68s
% Computational Cost: add. (299->91), mult. (664->155), div. (0->0), fcn. (793->12), ass. (0->43)
t113 = MDP(10) - MDP(13);
t112 = MDP(11) - MDP(14);
t79 = sin(pkin(11));
t83 = sin(qJ(2));
t86 = cos(qJ(2));
t96 = cos(pkin(11));
t97 = cos(pkin(6));
t92 = t97 * t96;
t65 = t79 * t83 - t86 * t92;
t95 = t79 * t97;
t67 = t96 * t83 + t86 * t95;
t111 = -g(1) * t67 - g(2) * t65;
t80 = sin(pkin(6));
t108 = g(3) * t80;
t78 = qJ(5) + qJ(6);
t76 = sin(t78);
t82 = sin(qJ(3));
t107 = t76 * t82;
t77 = cos(t78);
t106 = t77 * t82;
t105 = t80 * t83;
t85 = cos(qJ(3));
t104 = t80 * t85;
t103 = t80 * t86;
t81 = sin(qJ(5));
t102 = t81 * t82;
t84 = cos(qJ(5));
t101 = t82 * t84;
t100 = t82 * t86;
t99 = t84 * t86;
t66 = t79 * t86 + t83 * t92;
t94 = t80 * t96;
t61 = t66 * t82 + t85 * t94;
t68 = -t83 * t95 + t96 * t86;
t63 = -t79 * t104 + t68 * t82;
t69 = t82 * t105 - t97 * t85;
t98 = (-g(1) * (t63 * t77 - t67 * t76) - g(2) * (t61 * t77 - t65 * t76) - g(3) * (t76 * t103 + t69 * t77)) * MDP(28) + (-g(1) * (-t63 * t76 - t67 * t77) - g(2) * (-t61 * t76 - t65 * t77) - g(3) * (t77 * t103 - t69 * t76)) * MDP(29);
t93 = -g(1) * t68 - g(2) * t66;
t90 = g(1) * t63 + g(2) * t61 + g(3) * t69;
t70 = t83 * t104 + t97 * t82;
t64 = t79 * t80 * t82 + t68 * t85;
t62 = t66 * t85 - t82 * t94;
t1 = [(-MDP(1) - MDP(15)) * g(3); ((-t83 * t108 + t93) * pkin(8) + (-t86 * t108 - t111) * (pkin(3) * t85 + qJ(4) * t82 + pkin(2))) * MDP(15) + (-g(1) * (-t67 * t102 + t68 * t84) - g(2) * (-t65 * t102 + t66 * t84) - (t81 * t100 + t83 * t84) * t108) * MDP(21) + (-g(1) * (-t67 * t101 - t68 * t81) - g(2) * (-t65 * t101 - t66 * t81) - (-t81 * t83 + t82 * t99) * t108) * MDP(22) + (-g(1) * (-t67 * t107 + t68 * t77) - g(2) * (-t65 * t107 + t66 * t77) - (t76 * t100 + t77 * t83) * t108) * MDP(28) + (-g(1) * (-t67 * t106 - t68 * t76) - g(2) * (-t65 * t106 - t66 * t76) - (t77 * t100 - t76 * t83) * t108) * MDP(29) + (MDP(4) - MDP(12)) * (g(3) * t105 - t93) + (t112 * t82 - t113 * t85 - MDP(3)) * (g(3) * t103 + t111); (-g(1) * (-pkin(3) * t63 + qJ(4) * t64) - g(2) * (-pkin(3) * t61 + qJ(4) * t62) - g(3) * (-pkin(3) * t69 + qJ(4) * t70)) * MDP(15) + (-t81 * MDP(21) - t84 * MDP(22) - t76 * MDP(28) - t77 * MDP(29) + t112) * (g(1) * t64 + g(2) * t62 + g(3) * t70) + t113 * t90; -t90 * MDP(15); (-g(1) * (t63 * t84 - t67 * t81) - g(2) * (t61 * t84 - t65 * t81) - g(3) * (t81 * t103 + t69 * t84)) * MDP(21) + (-g(1) * (-t63 * t81 - t67 * t84) - g(2) * (-t61 * t81 - t65 * t84) - g(3) * (-t69 * t81 + t80 * t99)) * MDP(22) + t98; t98;];
taug  = t1;
