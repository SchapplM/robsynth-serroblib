% Calculate Gravitation load on the joints for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:38:57
% EndTime: 2019-03-08 20:38:58
% DurationCPUTime: 0.36s
% Computational Cost: add. (344->85), mult. (552->150), div. (0->0), fcn. (660->14), ass. (0->38)
t70 = sin(pkin(6));
t95 = g(3) * t70;
t66 = pkin(12) + qJ(4);
t63 = cos(t66);
t67 = qJ(5) + qJ(6);
t64 = sin(t67);
t94 = t63 * t64;
t65 = cos(t67);
t93 = t63 * t65;
t72 = sin(qJ(5));
t92 = t63 * t72;
t74 = cos(qJ(5));
t91 = t63 * t74;
t75 = cos(qJ(2));
t90 = t63 * t75;
t69 = sin(pkin(11));
t89 = t69 * t70;
t73 = sin(qJ(2));
t88 = t70 * t73;
t87 = t70 * t75;
t86 = t72 * t75;
t85 = t74 * t75;
t82 = cos(pkin(11));
t83 = cos(pkin(6));
t79 = t83 * t82;
t57 = t69 * t75 + t73 * t79;
t62 = sin(t66);
t80 = t70 * t82;
t51 = t57 * t63 - t62 * t80;
t81 = t69 * t83;
t59 = -t73 * t81 + t82 * t75;
t53 = t59 * t63 + t62 * t89;
t55 = t83 * t62 + t63 * t88;
t56 = t69 * t73 - t75 * t79;
t58 = t82 * t73 + t75 * t81;
t84 = (-g(1) * (-t53 * t64 + t58 * t65) - g(2) * (-t51 * t64 + t56 * t65) - g(3) * (-t55 * t64 - t65 * t87)) * MDP(28) + (-g(1) * (-t53 * t65 - t58 * t64) - g(2) * (-t51 * t65 - t56 * t64) - g(3) * (-t55 * t65 + t64 * t87)) * MDP(29);
t77 = -g(1) * t58 - g(2) * t56 + g(3) * t87;
t1 = [(-MDP(1) - MDP(8)) * g(3); (-g(1) * (-t58 * pkin(2) + t59 * qJ(3)) - g(2) * (-t56 * pkin(2) + t57 * qJ(3)) - (pkin(2) * t75 + qJ(3) * t73) * t95) * MDP(8) + (-g(1) * (-t58 * t91 + t59 * t72) - g(2) * (-t56 * t91 + t57 * t72) - (t63 * t85 + t72 * t73) * t95) * MDP(21) + (-g(1) * (t58 * t92 + t59 * t74) - g(2) * (t56 * t92 + t57 * t74) - (-t63 * t86 + t73 * t74) * t95) * MDP(22) + (-g(1) * (-t58 * t93 + t59 * t64) - g(2) * (-t56 * t93 + t57 * t64) - (t64 * t73 + t65 * t90) * t95) * MDP(28) + (-g(1) * (t58 * t94 + t59 * t65) - g(2) * (t56 * t94 + t57 * t65) - (-t64 * t90 + t65 * t73) * t95) * MDP(29) + (MDP(4) - MDP(7)) * (g(1) * t59 + g(2) * t57 + g(3) * t88) + (-t63 * MDP(14) + t62 * MDP(15) - MDP(5) * cos(pkin(12)) + MDP(6) * sin(pkin(12)) - MDP(3)) * t77; t77 * MDP(8); (g(1) * t53 + g(2) * t51 + g(3) * t55) * MDP(15) + (-MDP(21) * t74 + MDP(22) * t72 - MDP(28) * t65 + MDP(29) * t64 - MDP(14)) * (g(1) * (-t59 * t62 + t63 * t89) + g(2) * (-t57 * t62 - t63 * t80) + g(3) * (-t62 * t88 + t83 * t63)); (-g(1) * (-t53 * t72 + t58 * t74) - g(2) * (-t51 * t72 + t56 * t74) - g(3) * (-t55 * t72 - t70 * t85)) * MDP(21) + (-g(1) * (-t53 * t74 - t58 * t72) - g(2) * (-t51 * t74 - t56 * t72) - g(3) * (-t55 * t74 + t70 * t86)) * MDP(22) + t84; t84;];
taug  = t1;
