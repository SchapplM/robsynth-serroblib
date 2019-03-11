% Calculate Gravitation load on the joints for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:20:29
% EndTime: 2019-03-08 23:20:30
% DurationCPUTime: 0.50s
% Computational Cost: add. (353->95), mult. (677->160), div. (0->0), fcn. (805->12), ass. (0->45)
t108 = MDP(11) - MDP(19);
t80 = sin(qJ(2));
t83 = cos(qJ(2));
t95 = cos(pkin(11));
t96 = cos(pkin(6));
t90 = t96 * t95;
t94 = sin(pkin(11));
t61 = t94 * t80 - t83 * t90;
t107 = g(2) * t61;
t62 = t80 * t90 + t94 * t83;
t106 = g(2) * t62;
t76 = sin(pkin(6));
t105 = g(3) * t76;
t75 = qJ(4) + pkin(12) + qJ(6);
t72 = sin(t75);
t82 = cos(qJ(3));
t104 = t72 * t82;
t73 = cos(t75);
t103 = t73 * t82;
t102 = t76 * t80;
t101 = t76 * t83;
t78 = sin(qJ(4));
t100 = t78 * t82;
t81 = cos(qJ(4));
t99 = t81 * t82;
t98 = t82 * t83;
t79 = sin(qJ(3));
t92 = t76 * t95;
t58 = t62 * t82 - t79 * t92;
t89 = t96 * t94;
t64 = -t80 * t89 + t95 * t83;
t91 = t76 * t94;
t60 = t64 * t82 + t79 * t91;
t63 = t95 * t80 + t83 * t89;
t66 = t82 * t102 + t96 * t79;
t97 = (-g(1) * (-t60 * t72 + t63 * t73) - g(2) * (-t58 * t72 + t61 * t73) - g(3) * (-t73 * t101 - t66 * t72)) * MDP(26) + (-g(1) * (-t60 * t73 - t63 * t72) - g(2) * (-t58 * t73 - t61 * t72) - g(3) * (t72 * t101 - t66 * t73)) * MDP(27);
t93 = pkin(4) * t78 + pkin(8);
t74 = pkin(4) * t81 + pkin(3);
t77 = -qJ(5) - pkin(9);
t88 = t74 * t82 - t77 * t79 + pkin(2);
t57 = t62 * t79 + t82 * t92;
t59 = t64 * t79 - t82 * t91;
t65 = t79 * t102 - t96 * t82;
t87 = g(1) * t59 + g(2) * t57 + g(3) * t65;
t1 = [(-MDP(1) - MDP(20)) * g(3); (g(1) * t64 + g(3) * t102 + t106) * MDP(4) + (-g(1) * (-t63 * t99 + t64 * t78) - g(2) * (-t61 * t99 + t62 * t78) - (t78 * t80 + t81 * t98) * t105) * MDP(17) + (-g(1) * (t63 * t100 + t64 * t81) - g(2) * (t61 * t100 + t62 * t81) - (-t78 * t98 + t80 * t81) * t105) * MDP(18) + (-g(1) * (-t88 * t63 + t93 * t64) - t93 * t106 + t88 * t107 - (t93 * t80 + t88 * t83) * t105) * MDP(20) + (-g(1) * (-t63 * t103 + t64 * t72) - g(2) * (-t61 * t103 + t62 * t72) - (t72 * t80 + t73 * t98) * t105) * MDP(26) + (-g(1) * (t63 * t104 + t64 * t73) - g(2) * (t61 * t104 + t62 * t73) - (-t72 * t98 + t73 * t80) * t105) * MDP(27) + (-t82 * MDP(10) + t108 * t79 - MDP(3)) * (-g(1) * t63 + g(3) * t101 - t107); (-g(1) * (-t59 * t74 - t60 * t77) - g(2) * (-t57 * t74 - t58 * t77) - g(3) * (-t65 * t74 - t66 * t77)) * MDP(20) + t108 * (g(1) * t60 + g(2) * t58 + g(3) * t66) + (MDP(17) * t81 - MDP(18) * t78 + MDP(26) * t73 - MDP(27) * t72 + MDP(10)) * t87; (-g(1) * (-t60 * t81 - t63 * t78) - g(2) * (-t58 * t81 - t61 * t78) - g(3) * (t78 * t101 - t66 * t81)) * MDP(18) + t97 + (pkin(4) * MDP(20) + MDP(17)) * (-g(1) * (-t60 * t78 + t63 * t81) - g(2) * (-t58 * t78 + t61 * t81) - g(3) * (-t81 * t101 - t66 * t78)); -t87 * MDP(20); t97;];
taug  = t1;
