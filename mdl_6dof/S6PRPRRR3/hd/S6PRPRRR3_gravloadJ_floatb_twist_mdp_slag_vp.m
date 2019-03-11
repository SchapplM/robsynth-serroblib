% Calculate Gravitation load on the joints for
% S6PRPRRR3
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
%   see S6PRPRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:34:13
% EndTime: 2019-03-08 20:34:14
% DurationCPUTime: 0.34s
% Computational Cost: add. (344->73), mult. (474->126), div. (0->0), fcn. (552->14), ass. (0->34)
t70 = sin(pkin(6));
t91 = g(3) * t70;
t67 = pkin(12) + qJ(4);
t66 = qJ(5) + t67;
t63 = cos(t66);
t72 = sin(qJ(6));
t90 = t63 * t72;
t74 = cos(qJ(6));
t89 = t63 * t74;
t69 = sin(pkin(11));
t88 = t69 * t70;
t73 = sin(qJ(2));
t87 = t70 * t73;
t75 = cos(qJ(2));
t86 = t72 * t75;
t85 = t74 * t75;
t84 = cos(pkin(6));
t83 = cos(pkin(11));
t82 = t69 * t84;
t81 = t70 * t83;
t79 = t84 * t83;
t57 = t69 * t75 + t73 * t79;
t62 = sin(t66);
t51 = t57 * t63 - t62 * t81;
t59 = -t73 * t82 + t83 * t75;
t53 = t59 * t63 + t62 * t88;
t55 = t84 * t62 + t63 * t87;
t80 = (g(1) * t53 + g(2) * t51 + g(3) * t55) * MDP(22) + (-t74 * MDP(28) + t72 * MDP(29) - MDP(21)) * (g(1) * (-t59 * t62 + t63 * t88) + g(2) * (-t57 * t62 - t63 * t81) + g(3) * (-t62 * t87 + t84 * t63));
t56 = t69 * t73 - t75 * t79;
t58 = t83 * t73 + t75 * t82;
t77 = -g(1) * t58 - g(2) * t56 + t75 * t91;
t65 = cos(t67);
t64 = sin(t67);
t1 = [(-MDP(1) - MDP(8)) * g(3); (-g(1) * (-t58 * pkin(2) + t59 * qJ(3)) - g(2) * (-t56 * pkin(2) + t57 * qJ(3)) - (pkin(2) * t75 + qJ(3) * t73) * t91) * MDP(8) + (-g(1) * (-t58 * t89 + t59 * t72) - g(2) * (-t56 * t89 + t57 * t72) - (t63 * t85 + t72 * t73) * t91) * MDP(28) + (-g(1) * (t58 * t90 + t59 * t74) - g(2) * (t56 * t90 + t57 * t74) - (-t63 * t86 + t73 * t74) * t91) * MDP(29) + (MDP(4) - MDP(7)) * (g(1) * t59 + g(2) * t57 + g(3) * t87) + (-MDP(14) * t65 + MDP(15) * t64 - t63 * MDP(21) + MDP(22) * t62 - MDP(5) * cos(pkin(12)) + MDP(6) * sin(pkin(12)) - MDP(3)) * t77; t77 * MDP(8); (-g(1) * (-t59 * t64 + t65 * t88) - g(2) * (-t57 * t64 - t65 * t81) - g(3) * (-t64 * t87 + t84 * t65)) * MDP(14) + (-g(1) * (-t59 * t65 - t64 * t88) - g(2) * (-t57 * t65 + t64 * t81) - g(3) * (-t84 * t64 - t65 * t87)) * MDP(15) + t80; t80; (-g(1) * (-t53 * t72 + t58 * t74) - g(2) * (-t51 * t72 + t56 * t74) - g(3) * (-t55 * t72 - t70 * t85)) * MDP(28) + (-g(1) * (-t53 * t74 - t58 * t72) - g(2) * (-t51 * t74 - t56 * t72) - g(3) * (-t55 * t74 + t70 * t86)) * MDP(29);];
taug  = t1;
