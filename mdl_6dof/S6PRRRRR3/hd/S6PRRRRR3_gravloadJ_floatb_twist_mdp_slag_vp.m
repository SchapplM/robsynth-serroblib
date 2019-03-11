% Calculate Gravitation load on the joints for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:51:33
% EndTime: 2019-03-09 00:51:34
% DurationCPUTime: 0.44s
% Computational Cost: add. (446->102), mult. (767->180), div. (0->0), fcn. (946->14), ass. (0->40)
t70 = sin(pkin(6));
t96 = g(3) * t70;
t68 = qJ(4) + qJ(5);
t67 = qJ(6) + t68;
t63 = sin(t67);
t75 = cos(qJ(3));
t95 = t63 * t75;
t64 = cos(t67);
t94 = t64 * t75;
t65 = sin(t68);
t93 = t65 * t75;
t66 = cos(t68);
t92 = t66 * t75;
t73 = sin(qJ(2));
t91 = t70 * t73;
t90 = t70 * t75;
t76 = cos(qJ(2));
t89 = t70 * t76;
t71 = sin(qJ(4));
t88 = t71 * t75;
t74 = cos(qJ(4));
t87 = t74 * t75;
t86 = t75 * t76;
t69 = sin(pkin(12));
t83 = cos(pkin(12));
t84 = cos(pkin(6));
t79 = t84 * t83;
t57 = t69 * t76 + t73 * t79;
t72 = sin(qJ(3));
t81 = t70 * t83;
t53 = t57 * t75 - t72 * t81;
t82 = t69 * t84;
t59 = -t73 * t82 + t76 * t83;
t55 = t69 * t70 * t72 + t59 * t75;
t56 = t69 * t73 - t76 * t79;
t58 = t73 * t83 + t76 * t82;
t61 = t72 * t84 + t73 * t90;
t85 = (-g(1) * (-t55 * t63 + t58 * t64) - g(2) * (-t53 * t63 + t56 * t64) - g(3) * (-t61 * t63 - t64 * t89)) * MDP(31) + (-g(1) * (-t55 * t64 - t58 * t63) - g(2) * (-t53 * t64 - t56 * t63) - g(3) * (-t61 * t64 + t63 * t89)) * MDP(32);
t80 = (-g(1) * (-t55 * t65 + t58 * t66) - g(2) * (-t53 * t65 + t56 * t66) - g(3) * (-t61 * t65 - t66 * t89)) * MDP(24) + (-g(1) * (-t55 * t66 - t58 * t65) - g(2) * (-t53 * t66 - t56 * t65) - g(3) * (-t61 * t66 + t65 * t89)) * MDP(25) + t85;
t1 = [-g(3) * MDP(1); (g(1) * t59 + g(2) * t57 + g(3) * t91) * MDP(4) + (-g(1) * (-t58 * t87 + t59 * t71) - g(2) * (-t56 * t87 + t57 * t71) - (t71 * t73 + t74 * t86) * t96) * MDP(17) + (-g(1) * (t58 * t88 + t59 * t74) - g(2) * (t56 * t88 + t57 * t74) - (-t71 * t86 + t73 * t74) * t96) * MDP(18) + (-g(1) * (-t58 * t92 + t59 * t65) - g(2) * (-t56 * t92 + t57 * t65) - (t65 * t73 + t66 * t86) * t96) * MDP(24) + (-g(1) * (t58 * t93 + t59 * t66) - g(2) * (t56 * t93 + t57 * t66) - (-t65 * t86 + t66 * t73) * t96) * MDP(25) + (-g(1) * (-t58 * t94 + t59 * t63) - g(2) * (-t56 * t94 + t57 * t63) - (t63 * t73 + t64 * t86) * t96) * MDP(31) + (-g(1) * (t58 * t95 + t59 * t64) - g(2) * (t56 * t95 + t57 * t64) - (-t63 * t86 + t64 * t73) * t96) * MDP(32) + (-t75 * MDP(10) + t72 * MDP(11) - MDP(3)) * (-g(1) * t58 - g(2) * t56 + g(3) * t89); (g(1) * t55 + g(2) * t53 + g(3) * t61) * MDP(11) + (-MDP(17) * t74 + MDP(18) * t71 - MDP(24) * t66 + MDP(25) * t65 - MDP(31) * t64 + MDP(32) * t63 - MDP(10)) * (g(1) * (-t59 * t72 + t69 * t90) + g(2) * (-t57 * t72 - t75 * t81) + g(3) * (-t72 * t91 + t75 * t84)); (-g(1) * (-t55 * t71 + t58 * t74) - g(2) * (-t53 * t71 + t56 * t74) - g(3) * (-t61 * t71 - t74 * t89)) * MDP(17) + (-g(1) * (-t55 * t74 - t58 * t71) - g(2) * (-t53 * t74 - t56 * t71) - g(3) * (-t61 * t74 + t71 * t89)) * MDP(18) + t80; t80; t85;];
taug  = t1;
