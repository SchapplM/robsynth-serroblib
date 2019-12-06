% Calculate Gravitation load on the joints for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:56
% EndTime: 2019-12-05 17:44:58
% DurationCPUTime: 0.27s
% Computational Cost: add. (154->56), mult. (178->86), div. (0->0), fcn. (175->10), ass. (0->40)
t63 = sin(pkin(8));
t83 = g(1) * t63;
t67 = cos(qJ(1));
t82 = g(2) * t67;
t61 = pkin(9) + qJ(4);
t59 = qJ(5) + t61;
t55 = sin(t59);
t66 = sin(qJ(1));
t81 = t66 * t55;
t56 = cos(t59);
t80 = t66 * t56;
t57 = sin(t61);
t79 = t66 * t57;
t58 = cos(t61);
t78 = t66 * t58;
t62 = sin(pkin(9));
t77 = t66 * t62;
t64 = cos(pkin(9));
t76 = t66 * t64;
t75 = t67 * t55;
t74 = t67 * t56;
t73 = t67 * t57;
t72 = t67 * t58;
t71 = t67 * t62;
t70 = t67 * t64;
t65 = cos(pkin(8));
t44 = t65 * t81 + t74;
t45 = t65 * t80 - t75;
t46 = t65 * t75 - t80;
t47 = -t65 * t74 - t81;
t69 = (-g(2) * t44 + g(3) * t46 + t55 * t83) * MDP(24) + (-g(2) * t45 - g(3) * t47 + t56 * t83) * MDP(25);
t54 = g(3) * t66 + t82;
t53 = g(2) * t66 - g(3) * t67;
t68 = pkin(2) * t65 + qJ(3) * t63 + pkin(1);
t60 = t67 * qJ(2);
t51 = -t65 * t72 - t79;
t50 = t65 * t73 - t78;
t49 = t65 * t78 - t73;
t48 = t65 * t79 + t72;
t1 = [(-g(2) * (-t67 * pkin(1) - t66 * qJ(2)) - g(3) * (-t66 * pkin(1) + t60)) * MDP(7) + (-g(2) * (-t65 * t70 - t77) - g(3) * (-t65 * t76 + t71)) * MDP(8) + (-g(2) * (t65 * t71 - t76) - g(3) * (t65 * t77 + t70)) * MDP(9) + (-g(3) * t60 + t68 * t82 + (g(2) * qJ(2) + g(3) * t68) * t66) * MDP(11) + (-g(2) * t51 + g(3) * t49) * MDP(17) + (-g(2) * t50 - g(3) * t48) * MDP(18) + (-g(2) * t47 + g(3) * t45) * MDP(24) + (-g(2) * t46 - g(3) * t44) * MDP(25) + (-MDP(3) + MDP(6)) * t53 + (t65 * MDP(4) + MDP(2) + (-MDP(5) + MDP(10)) * t63) * t54; (-MDP(11) - MDP(7)) * t54; (g(1) * t65 + t53 * t63) * MDP(11); (-g(2) * t48 + g(3) * t50 + t57 * t83) * MDP(17) + (-g(2) * t49 - g(3) * t51 + t58 * t83) * MDP(18) + t69; t69;];
taug = t1;
