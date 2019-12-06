% Calculate Gravitation load on the joints for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5PRRRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:13
% EndTime: 2019-12-05 17:10:14
% DurationCPUTime: 0.12s
% Computational Cost: add. (155->34), mult. (175->51), div. (0->0), fcn. (166->10), ass. (0->25)
t46 = qJ(2) + qJ(3);
t42 = sin(t46);
t67 = g(3) * t42;
t45 = qJ(4) + qJ(5);
t41 = sin(t45);
t47 = sin(pkin(9));
t65 = t47 * t41;
t43 = cos(t45);
t64 = t47 * t43;
t49 = sin(qJ(4));
t63 = t47 * t49;
t51 = cos(qJ(4));
t62 = t47 * t51;
t48 = cos(pkin(9));
t61 = t48 * t41;
t60 = t48 * t43;
t59 = t48 * t49;
t58 = t48 * t51;
t44 = cos(t46);
t57 = (-g(1) * (-t44 * t61 + t64) - g(2) * (-t44 * t65 - t60) + t41 * t67) * MDP(20) + (-g(1) * (-t44 * t60 - t65) - g(2) * (-t44 * t64 + t61) + t43 * t67) * MDP(21);
t56 = g(1) * t48 + g(2) * t47;
t55 = (t56 * t44 + t67) * MDP(7) + (t51 * MDP(13) - t49 * MDP(14) + t43 * MDP(20) - t41 * MDP(21) + MDP(6)) * (-g(3) * t44 + t56 * t42);
t52 = cos(qJ(2));
t50 = sin(qJ(2));
t1 = [-g(3) * MDP(1); (-g(3) * t52 + t56 * t50) * MDP(3) + (g(3) * t50 + t56 * t52) * MDP(4) + t55; t55; (-g(1) * (-t44 * t59 + t62) - g(2) * (-t44 * t63 - t58) + t49 * t67) * MDP(13) + (-g(1) * (-t44 * t58 - t63) - g(2) * (-t44 * t62 + t59) + t51 * t67) * MDP(14) + t57; t57;];
taug = t1;
