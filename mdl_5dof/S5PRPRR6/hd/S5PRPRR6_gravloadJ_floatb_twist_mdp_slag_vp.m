% Calculate Gravitation load on the joints for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRPRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:57:54
% EndTime: 2019-12-05 15:57:55
% DurationCPUTime: 0.25s
% Computational Cost: add. (173->58), mult. (326->104), div. (0->0), fcn. (378->12), ass. (0->30)
t56 = sin(pkin(5));
t76 = g(3) * t56;
t53 = pkin(10) + qJ(4);
t52 = cos(t53);
t58 = sin(qJ(5));
t75 = t52 * t58;
t60 = cos(qJ(5));
t74 = t52 * t60;
t55 = sin(pkin(9));
t73 = t55 * t56;
t59 = sin(qJ(2));
t72 = t56 * t59;
t61 = cos(qJ(2));
t71 = t58 * t61;
t70 = t60 * t61;
t69 = cos(pkin(5));
t68 = cos(pkin(9));
t67 = t55 * t69;
t66 = t56 * t68;
t65 = t69 * t68;
t45 = t55 * t59 - t61 * t65;
t47 = t68 * t59 + t61 * t67;
t63 = -g(1) * t47 - g(2) * t45 + t61 * t76;
t51 = sin(t53);
t48 = -t59 * t67 + t68 * t61;
t46 = t55 * t61 + t59 * t65;
t44 = t69 * t51 + t52 * t72;
t42 = t48 * t52 + t51 * t73;
t40 = t46 * t52 - t51 * t66;
t1 = [(-MDP(1) - MDP(8)) * g(3); (-g(1) * (-pkin(2) * t47 + qJ(3) * t48) - g(2) * (-pkin(2) * t45 + qJ(3) * t46) - (pkin(2) * t61 + qJ(3) * t59) * t76) * MDP(8) + (-g(1) * (-t47 * t74 + t48 * t58) - g(2) * (-t45 * t74 + t46 * t58) - (t52 * t70 + t58 * t59) * t76) * MDP(21) + (-g(1) * (t47 * t75 + t48 * t60) - g(2) * (t45 * t75 + t46 * t60) - (-t52 * t71 + t59 * t60) * t76) * MDP(22) + (MDP(4) - MDP(7)) * (g(1) * t48 + g(2) * t46 + g(3) * t72) + (-t52 * MDP(14) + MDP(15) * t51 - MDP(5) * cos(pkin(10)) + MDP(6) * sin(pkin(10)) - MDP(3)) * t63; t63 * MDP(8); (g(1) * t42 + g(2) * t40 + g(3) * t44) * MDP(15) + (-MDP(21) * t60 + MDP(22) * t58 - MDP(14)) * (g(1) * (-t48 * t51 + t52 * t73) + g(2) * (-t46 * t51 - t52 * t66) + g(3) * (-t51 * t72 + t69 * t52)); (-g(1) * (-t42 * t58 + t47 * t60) - g(2) * (-t40 * t58 + t45 * t60) - g(3) * (-t44 * t58 - t56 * t70)) * MDP(21) + (-g(1) * (-t42 * t60 - t47 * t58) - g(2) * (-t40 * t60 - t45 * t58) - g(3) * (-t44 * t60 + t56 * t71)) * MDP(22);];
taug = t1;
