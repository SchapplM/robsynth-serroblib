% Calculate Gravitation load on the joints for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:16:42
% EndTime: 2019-12-05 17:16:43
% DurationCPUTime: 0.29s
% Computational Cost: add. (225->63), mult. (399->115), div. (0->0), fcn. (474->12), ass. (0->34)
t62 = sin(pkin(5));
t85 = g(3) * t62;
t60 = qJ(3) + qJ(4);
t59 = cos(t60);
t65 = sin(qJ(5));
t84 = t59 * t65;
t68 = cos(qJ(5));
t83 = t59 * t68;
t61 = sin(pkin(10));
t82 = t61 * t62;
t63 = cos(pkin(10));
t81 = t62 * t63;
t66 = sin(qJ(3));
t80 = t62 * t66;
t67 = sin(qJ(2));
t79 = t62 * t67;
t69 = cos(qJ(3));
t78 = t62 * t69;
t64 = cos(pkin(5));
t77 = t64 * t67;
t70 = cos(qJ(2));
t76 = t64 * t70;
t75 = t65 * t70;
t74 = t68 * t70;
t54 = t61 * t70 + t63 * t77;
t58 = sin(t60);
t48 = t54 * t59 - t58 * t81;
t56 = -t61 * t77 + t63 * t70;
t50 = t56 * t59 + t58 * t82;
t52 = t64 * t58 + t59 * t79;
t73 = (g(1) * t50 + g(2) * t48 + g(3) * t52) * MDP(18) + (-MDP(24) * t68 + MDP(25) * t65 - MDP(17)) * (g(1) * (-t56 * t58 + t59 * t82) + g(2) * (-t54 * t58 - t59 * t81) + g(3) * (-t58 * t79 + t64 * t59));
t55 = t61 * t76 + t63 * t67;
t53 = t61 * t67 - t63 * t76;
t1 = [-g(3) * MDP(1); (g(1) * t56 + g(2) * t54 + g(3) * t79) * MDP(4) + (-g(1) * (-t55 * t83 + t56 * t65) - g(2) * (-t53 * t83 + t54 * t65) - (t59 * t74 + t65 * t67) * t85) * MDP(24) + (-g(1) * (t55 * t84 + t56 * t68) - g(2) * (t53 * t84 + t54 * t68) - (-t59 * t75 + t67 * t68) * t85) * MDP(25) + (-MDP(10) * t69 + MDP(11) * t66 - t59 * MDP(17) + MDP(18) * t58 - MDP(3)) * (-g(1) * t55 - g(2) * t53 + t70 * t85); (-g(1) * (-t56 * t66 + t61 * t78) - g(2) * (-t54 * t66 - t63 * t78) - g(3) * (t64 * t69 - t66 * t79)) * MDP(10) + (-g(1) * (-t56 * t69 - t61 * t80) - g(2) * (-t54 * t69 + t63 * t80) - g(3) * (-t64 * t66 - t67 * t78)) * MDP(11) + t73; t73; (-g(1) * (-t50 * t65 + t55 * t68) - g(2) * (-t48 * t65 + t53 * t68) - g(3) * (-t52 * t65 - t62 * t74)) * MDP(24) + (-g(1) * (-t50 * t68 - t55 * t65) - g(2) * (-t48 * t68 - t53 * t65) - g(3) * (-t52 * t68 + t62 * t75)) * MDP(25);];
taug = t1;
