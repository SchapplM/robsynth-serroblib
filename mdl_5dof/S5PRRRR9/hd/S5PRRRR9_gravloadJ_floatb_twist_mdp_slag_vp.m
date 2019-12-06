% Calculate Gravitation load on the joints for
% S5PRRRR9
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
%   see S5PRRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:21:11
% EndTime: 2019-12-05 17:21:13
% DurationCPUTime: 0.28s
% Computational Cost: add. (225->75), mult. (477->136), div. (0->0), fcn. (582->12), ass. (0->34)
t58 = sin(pkin(5));
t81 = g(3) * t58;
t56 = qJ(4) + qJ(5);
t54 = sin(t56);
t63 = cos(qJ(3));
t80 = t54 * t63;
t55 = cos(t56);
t79 = t55 * t63;
t61 = sin(qJ(2));
t78 = t58 * t61;
t77 = t58 * t63;
t64 = cos(qJ(2));
t76 = t58 * t64;
t59 = sin(qJ(4));
t75 = t59 * t63;
t62 = cos(qJ(4));
t74 = t62 * t63;
t73 = t63 * t64;
t57 = sin(pkin(10));
t70 = cos(pkin(10));
t71 = cos(pkin(5));
t67 = t71 * t70;
t48 = t57 * t64 + t61 * t67;
t60 = sin(qJ(3));
t68 = t58 * t70;
t44 = t48 * t63 - t60 * t68;
t69 = t57 * t71;
t50 = -t61 * t69 + t70 * t64;
t46 = t57 * t58 * t60 + t50 * t63;
t47 = t57 * t61 - t64 * t67;
t49 = t70 * t61 + t64 * t69;
t52 = t71 * t60 + t61 * t77;
t72 = (-g(1) * (-t46 * t54 + t49 * t55) - g(2) * (-t44 * t54 + t47 * t55) - g(3) * (-t52 * t54 - t55 * t76)) * MDP(24) + (-g(1) * (-t46 * t55 - t49 * t54) - g(2) * (-t44 * t55 - t47 * t54) - g(3) * (-t52 * t55 + t54 * t76)) * MDP(25);
t1 = [-g(3) * MDP(1); (g(1) * t50 + g(2) * t48 + g(3) * t78) * MDP(4) + (-g(1) * (-t49 * t74 + t50 * t59) - g(2) * (-t47 * t74 + t48 * t59) - (t59 * t61 + t62 * t73) * t81) * MDP(17) + (-g(1) * (t49 * t75 + t50 * t62) - g(2) * (t47 * t75 + t48 * t62) - (-t59 * t73 + t61 * t62) * t81) * MDP(18) + (-g(1) * (-t49 * t79 + t50 * t54) - g(2) * (-t47 * t79 + t48 * t54) - (t54 * t61 + t55 * t73) * t81) * MDP(24) + (-g(1) * (t49 * t80 + t50 * t55) - g(2) * (t47 * t80 + t48 * t55) - (-t54 * t73 + t55 * t61) * t81) * MDP(25) + (-t63 * MDP(10) + t60 * MDP(11) - MDP(3)) * (-g(1) * t49 - g(2) * t47 + g(3) * t76); (g(1) * t46 + g(2) * t44 + g(3) * t52) * MDP(11) + (-MDP(17) * t62 + MDP(18) * t59 - MDP(24) * t55 + MDP(25) * t54 - MDP(10)) * (g(1) * (-t50 * t60 + t57 * t77) + g(2) * (-t48 * t60 - t63 * t68) + g(3) * (-t60 * t78 + t71 * t63)); (-g(1) * (-t46 * t59 + t49 * t62) - g(2) * (-t44 * t59 + t47 * t62) - g(3) * (-t52 * t59 - t62 * t76)) * MDP(17) + (-g(1) * (-t46 * t62 - t49 * t59) - g(2) * (-t44 * t62 - t47 * t59) - g(3) * (-t52 * t62 + t59 * t76)) * MDP(18) + t72; t72;];
taug = t1;
