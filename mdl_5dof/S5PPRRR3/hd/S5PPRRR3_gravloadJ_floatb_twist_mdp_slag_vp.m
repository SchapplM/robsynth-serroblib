% Calculate Gravitation load on the joints for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PPRRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:17:00
% EndTime: 2019-12-05 15:17:01
% DurationCPUTime: 0.14s
% Computational Cost: add. (107->39), mult. (197->70), div. (0->0), fcn. (224->10), ass. (0->24)
t37 = sin(pkin(9));
t38 = sin(pkin(8));
t55 = t37 * t38;
t40 = cos(pkin(8));
t54 = t37 * t40;
t41 = sin(qJ(4));
t53 = t37 * t41;
t43 = cos(qJ(4));
t52 = t37 * t43;
t44 = cos(qJ(3));
t51 = t37 * t44;
t42 = sin(qJ(3));
t50 = t38 * t42;
t49 = t38 * t44;
t48 = t40 * t42;
t47 = t40 * t44;
t39 = cos(pkin(9));
t31 = t39 * t49 - t48;
t33 = t39 * t47 + t50;
t36 = qJ(4) + qJ(5);
t34 = sin(t36);
t35 = cos(t36);
t46 = (-g(1) * (-t33 * t34 + t35 * t54) - g(2) * (-t31 * t34 + t35 * t55) - g(3) * (-t34 * t51 - t35 * t39)) * MDP(18) + (-g(1) * (-t33 * t35 - t34 * t54) - g(2) * (-t31 * t35 - t34 * t55) - g(3) * (t34 * t39 - t35 * t51)) * MDP(19);
t1 = [(-MDP(1) - MDP(2)) * g(3); (-g(1) * t38 + g(2) * t40) * MDP(2); (g(1) * t33 + g(2) * t31 + g(3) * t51) * MDP(5) + (MDP(11) * t43 - MDP(12) * t41 + MDP(18) * t35 - MDP(19) * t34 + MDP(4)) * (g(3) * t37 * t42 - g(1) * (-t39 * t48 + t49) - g(2) * (-t39 * t50 - t47)); (-g(1) * (-t33 * t41 + t40 * t52) - g(2) * (-t31 * t41 + t38 * t52) - g(3) * (-t39 * t43 - t41 * t51)) * MDP(11) + (-g(1) * (-t33 * t43 - t40 * t53) - g(2) * (-t31 * t43 - t38 * t53) - g(3) * (t39 * t41 - t43 * t51)) * MDP(12) + t46; t46;];
taug = t1;
