% Calculate Gravitation load on the joints for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:49:35
% EndTime: 2019-03-08 19:49:37
% DurationCPUTime: 0.69s
% Computational Cost: add. (253->90), mult. (571->148), div. (0->0), fcn. (663->12), ass. (0->45)
t112 = MDP(14) - MDP(17);
t81 = sin(pkin(10));
t85 = sin(qJ(2));
t87 = cos(qJ(2));
t98 = cos(pkin(10));
t99 = cos(pkin(6));
t92 = t99 * t98;
t67 = t81 * t87 + t85 * t92;
t95 = t81 * t99;
t69 = -t85 * t95 + t98 * t87;
t111 = -g(1) * t69 - g(2) * t67;
t82 = sin(pkin(6));
t108 = g(3) * t82;
t79 = pkin(11) + qJ(6);
t77 = sin(t79);
t84 = sin(qJ(4));
t107 = t77 * t84;
t78 = cos(t79);
t106 = t78 * t84;
t80 = sin(pkin(11));
t105 = t80 * t84;
t104 = t81 * t82;
t103 = t82 * t85;
t102 = t82 * t87;
t83 = cos(pkin(11));
t101 = t83 * t84;
t100 = t84 * t85;
t97 = MDP(18) + MDP(7);
t96 = g(3) * (pkin(2) * t102 + qJ(3) * t103);
t94 = t82 * t98;
t86 = cos(qJ(4));
t93 = pkin(4) * t84 - qJ(5) * t86;
t68 = t98 * t85 + t87 * t95;
t58 = t84 * t104 - t68 * t86;
t66 = t81 * t85 - t87 * t92;
t60 = t66 * t86 + t84 * t94;
t70 = t86 * t102 + t99 * t84;
t90 = g(1) * t58 - g(2) * t60 + g(3) * t70;
t57 = -g(1) * t68 - g(2) * t66 + g(3) * t102;
t71 = -t84 * t102 + t99 * t86;
t65 = t68 * pkin(2);
t64 = t66 * pkin(2);
t61 = -t66 * t84 + t86 * t94;
t59 = t86 * t104 + t68 * t84;
t1 = [(-MDP(1) - t97) * g(3); (-g(1) * (qJ(3) * t69 - t65) - g(2) * (qJ(3) * t67 - t64) - t96) * MDP(7) + (-g(1) * (t69 * t101 - t68 * t80) - g(2) * (t67 * t101 - t66 * t80) - (t83 * t100 + t80 * t87) * t108) * MDP(15) + (-g(1) * (-t69 * t105 - t68 * t83) - g(2) * (-t67 * t105 - t66 * t83) - (-t80 * t100 + t83 * t87) * t108) * MDP(16) + (-g(1) * (-pkin(8) * t68 - t65) - g(2) * (-pkin(8) * t66 - t64) - t96 - (pkin(8) * t87 + t93 * t85) * t108 + t111 * (qJ(3) + t93)) * MDP(18) + (-g(1) * (t69 * t106 - t68 * t77) - g(2) * (t67 * t106 - t66 * t77) - (t78 * t100 + t77 * t87) * t108) * MDP(24) + (-g(1) * (-t69 * t107 - t68 * t78) - g(2) * (-t67 * t107 - t66 * t78) - (-t77 * t100 + t78 * t87) * t108) * MDP(25) + (-MDP(3) + MDP(5)) * t57 + (-t84 * MDP(13) - t112 * t86 + MDP(4) - MDP(6)) * (g(3) * t103 - t111); t97 * t57; (-g(1) * (-pkin(4) * t58 + qJ(5) * t59) - g(2) * (pkin(4) * t60 - qJ(5) * t61) - g(3) * (-pkin(4) * t70 + qJ(5) * t71)) * MDP(18) + t112 * (g(1) * t59 - g(2) * t61 + g(3) * t71) + (t83 * MDP(15) - t80 * MDP(16) + t78 * MDP(24) - t77 * MDP(25) + MDP(13)) * t90; -t90 * MDP(18); (-g(1) * (-t59 * t77 + t69 * t78) - g(2) * (t61 * t77 + t67 * t78) - g(3) * (t78 * t103 - t71 * t77)) * MDP(24) + (-g(1) * (-t59 * t78 - t69 * t77) - g(2) * (t61 * t78 - t67 * t77) - g(3) * (-t77 * t103 - t71 * t78)) * MDP(25);];
taug  = t1;
