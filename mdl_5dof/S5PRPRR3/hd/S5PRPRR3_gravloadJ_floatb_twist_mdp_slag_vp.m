% Calculate Gravitation load on the joints for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:47:33
% EndTime: 2019-12-05 15:47:34
% DurationCPUTime: 0.14s
% Computational Cost: add. (106->33), mult. (128->52), div. (0->0), fcn. (121->10), ass. (0->25)
t60 = MDP(5) * pkin(2) + MDP(3);
t36 = qJ(4) + qJ(5);
t33 = sin(t36);
t34 = cos(t36);
t39 = sin(qJ(4));
t41 = cos(qJ(4));
t59 = -MDP(11) * t41 + MDP(12) * t39 - MDP(18) * t34 + MDP(19) * t33;
t35 = qJ(2) + pkin(9);
t31 = sin(t35);
t58 = g(3) * t31;
t37 = sin(pkin(8));
t56 = t37 * t33;
t55 = t37 * t34;
t54 = t37 * t39;
t53 = t37 * t41;
t38 = cos(pkin(8));
t52 = t38 * t33;
t51 = t38 * t34;
t50 = t38 * t39;
t49 = t38 * t41;
t32 = cos(t35);
t48 = (-g(1) * (-t32 * t52 + t55) - g(2) * (-t32 * t56 - t51) + t33 * t58) * MDP(18) + (-g(1) * (-t32 * t51 - t56) - g(2) * (-t32 * t55 + t52) + t34 * t58) * MDP(19);
t42 = cos(qJ(2));
t40 = sin(qJ(2));
t1 = [(-MDP(1) - MDP(5)) * g(3); (t40 * MDP(4) + t59 * t32 - t60 * t42) * g(3) + (MDP(4) * t42 - t59 * t31 + t60 * t40) * (g(1) * t38 + g(2) * t37); (-g(1) * t37 + g(2) * t38) * MDP(5); (-g(1) * (-t32 * t50 + t53) - g(2) * (-t32 * t54 - t49) + t39 * t58) * MDP(11) + (-g(1) * (-t32 * t49 - t54) - g(2) * (-t32 * t53 + t50) + t41 * t58) * MDP(12) + t48; t48;];
taug = t1;
