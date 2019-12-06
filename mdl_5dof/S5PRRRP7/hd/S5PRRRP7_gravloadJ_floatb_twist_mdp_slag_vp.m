% Calculate Gravitation load on the joints for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRRP7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:26
% EndTime: 2019-12-05 16:56:28
% DurationCPUTime: 0.35s
% Computational Cost: add. (182->68), mult. (451->117), div. (0->0), fcn. (523->10), ass. (0->39)
t93 = MDP(11) - MDP(19);
t64 = sin(pkin(9));
t69 = sin(qJ(2));
t72 = cos(qJ(2));
t82 = cos(pkin(9));
t83 = cos(pkin(5));
t78 = t83 * t82;
t52 = t64 * t69 - t72 * t78;
t92 = g(2) * t52;
t53 = t64 * t72 + t69 * t78;
t91 = g(2) * t53;
t65 = sin(pkin(5));
t90 = g(3) * t65;
t89 = t65 * t69;
t71 = cos(qJ(3));
t88 = t65 * t71;
t87 = t65 * t72;
t67 = sin(qJ(4));
t86 = t67 * t71;
t70 = cos(qJ(4));
t85 = t70 * t71;
t84 = t71 * t72;
t81 = pkin(4) * t67 + pkin(7);
t80 = t64 * t83;
t79 = t65 * t82;
t63 = t70 * pkin(4) + pkin(3);
t66 = -qJ(5) - pkin(8);
t68 = sin(qJ(3));
t77 = t63 * t71 - t66 * t68 + pkin(2);
t48 = t53 * t68 + t71 * t79;
t55 = -t69 * t80 + t82 * t72;
t50 = t55 * t68 - t64 * t88;
t56 = t68 * t89 - t83 * t71;
t76 = g(1) * t50 + g(2) * t48 + g(3) * t56;
t57 = t83 * t68 + t69 * t88;
t54 = t82 * t69 + t72 * t80;
t51 = t64 * t65 * t68 + t55 * t71;
t49 = t53 * t71 - t68 * t79;
t1 = [(-MDP(1) - MDP(20)) * g(3); (g(1) * t55 + g(3) * t89 + t91) * MDP(4) + (-g(1) * (-t54 * t85 + t55 * t67) - g(2) * (-t52 * t85 + t53 * t67) - (t67 * t69 + t70 * t84) * t90) * MDP(17) + (-g(1) * (t54 * t86 + t55 * t70) - g(2) * (t52 * t86 + t53 * t70) - (-t67 * t84 + t69 * t70) * t90) * MDP(18) + (-g(1) * (-t77 * t54 + t81 * t55) - t81 * t91 + t77 * t92 - (t81 * t69 + t77 * t72) * t90) * MDP(20) + (-t71 * MDP(10) + t93 * t68 - MDP(3)) * (-g(1) * t54 + g(3) * t87 - t92); (-g(1) * (-t50 * t63 - t51 * t66) - g(2) * (-t48 * t63 - t49 * t66) - g(3) * (-t56 * t63 - t57 * t66)) * MDP(20) + t93 * (g(1) * t51 + g(2) * t49 + g(3) * t57) + (t70 * MDP(17) - t67 * MDP(18) + MDP(10)) * t76; (-g(1) * (-t51 * t70 - t54 * t67) - g(2) * (-t49 * t70 - t52 * t67) - g(3) * (-t57 * t70 + t67 * t87)) * MDP(18) + (pkin(4) * MDP(20) + MDP(17)) * (-g(1) * (-t51 * t67 + t54 * t70) - g(2) * (-t49 * t67 + t52 * t70) - g(3) * (-t57 * t67 - t70 * t87)); -t76 * MDP(20);];
taug = t1;
