% Calculate Gravitation load on the joints for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:26:43
% EndTime: 2019-03-09 15:26:44
% DurationCPUTime: 0.29s
% Computational Cost: add. (285->73), mult. (282->99), div. (0->0), fcn. (237->10), ass. (0->43)
t81 = sin(qJ(1));
t84 = cos(qJ(1));
t66 = g(1) * t84 + g(2) * t81;
t78 = qJ(2) + qJ(3);
t73 = sin(t78);
t100 = t66 * t73;
t72 = pkin(10) + t78;
t69 = sin(t72);
t70 = cos(t72);
t87 = t70 * pkin(4) + t69 * qJ(5);
t99 = MDP(19) + MDP(23);
t98 = pkin(4) * t69;
t97 = g(3) * t70;
t74 = cos(t78);
t96 = g(3) * t74;
t79 = sin(qJ(6));
t95 = t81 * t79;
t82 = cos(qJ(6));
t94 = t81 * t82;
t93 = t84 * t79;
t92 = t84 * t82;
t71 = pkin(3) * t74;
t83 = cos(qJ(2));
t75 = t83 * pkin(2);
t91 = t71 + t75;
t90 = qJ(5) * t70;
t89 = t71 + t87;
t80 = sin(qJ(2));
t62 = -t80 * pkin(2) - pkin(3) * t73;
t88 = t62 - t98;
t65 = g(1) * t81 - g(2) * t84;
t53 = -t66 * t69 + t97;
t86 = t53 * MDP(21) + (-t96 + t100) * MDP(16) + (g(3) * t73 + t66 * t74) * MDP(17) + (t79 * MDP(29) + t82 * MDP(30) + MDP(22)) * (-g(3) * t69 - t66 * t70);
t77 = -qJ(4) - pkin(8) - pkin(7);
t64 = t84 * t90;
t63 = t81 * t90;
t61 = pkin(1) + t91;
t60 = t84 * t61;
t59 = -t69 * t95 + t92;
t58 = t69 * t94 + t93;
t57 = t69 * t93 + t94;
t56 = t69 * t92 - t95;
t1 = [(-g(1) * (-t81 * t61 - t84 * t77) - g(2) * (-t81 * t77 + t60)) * MDP(19) + (-g(2) * t60 + (g(1) * t77 - g(2) * t87) * t84 + (-g(1) * (-t61 - t87) + g(2) * t77) * t81) * MDP(23) + (-g(1) * t59 - g(2) * t57) * MDP(29) + (g(1) * t58 - g(2) * t56) * MDP(30) + (MDP(3) - MDP(18) - MDP(20)) * t66 + (-t80 * MDP(10) + MDP(16) * t74 - MDP(17) * t73 - MDP(21) * t70 + MDP(22) * t69 + t83 * MDP(9) + MDP(2)) * t65; (-g(3) * t83 + t66 * t80) * MDP(9) + (g(3) * t80 + t66 * t83) * MDP(10) + (-g(3) * t91 - t62 * t66) * MDP(19) + (-g(1) * (t84 * t88 + t64) - g(2) * (t81 * t88 + t63) - g(3) * (t75 + t89)) * MDP(23) + t86; (-g(1) * (-t84 * t98 + t64) - g(2) * (-t81 * t98 + t63) - g(3) * t89) * MDP(23) + (-MDP(19) * t96 + t99 * t100) * pkin(3) + t86; -t99 * t65; t53 * MDP(23); (-g(1) * t56 - g(2) * t58 + t82 * t97) * MDP(29) + (g(1) * t57 - g(2) * t59 - t79 * t97) * MDP(30);];
taug  = t1;
