% Calculate Gravitation load on the joints for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:10
% EndTime: 2019-12-05 15:22:11
% DurationCPUTime: 0.19s
% Computational Cost: add. (127->43), mult. (119->68), div. (0->0), fcn. (107->8), ass. (0->27)
t46 = pkin(7) + qJ(2);
t42 = sin(t46);
t59 = g(1) * t42;
t48 = sin(pkin(8));
t58 = g(3) * t48;
t50 = cos(pkin(8));
t57 = t42 * t50;
t44 = cos(t46);
t56 = t44 * t50;
t47 = sin(pkin(9));
t55 = t47 * t50;
t49 = cos(pkin(9));
t54 = t49 * t50;
t53 = t44 * pkin(2) + t42 * qJ(3);
t52 = -MDP(12) - MDP(8);
t36 = g(1) * t44 + g(2) * t42;
t35 = -g(2) * t44 + t59;
t51 = pkin(3) * t50 + qJ(4) * t48;
t45 = pkin(9) + qJ(5);
t43 = cos(t45);
t41 = sin(t45);
t38 = t44 * qJ(3);
t33 = t41 * t42 + t43 * t56;
t32 = -t41 * t56 + t42 * t43;
t31 = t41 * t44 - t43 * t57;
t30 = t41 * t57 + t43 * t44;
t1 = [(-MDP(1) + t52) * g(3); (-g(1) * (-pkin(2) * t42 + t38) - g(2) * t53) * MDP(8) + (-g(1) * (-t42 * t54 + t44 * t47) - g(2) * (t42 * t47 + t44 * t54)) * MDP(9) + (-g(1) * (t42 * t55 + t44 * t49) - g(2) * (t42 * t49 - t44 * t55)) * MDP(10) + (-g(1) * t38 - g(2) * (t44 * t51 + t53) - (-pkin(2) - t51) * t59) * MDP(12) + (-g(1) * t31 - g(2) * t33) * MDP(18) + (-g(1) * t30 - g(2) * t32) * MDP(19) + (MDP(4) - MDP(7)) * t36 + (t50 * MDP(5) + MDP(3) + (-MDP(6) + MDP(11)) * t48) * t35; t52 * t35; (g(3) * t50 - t36 * t48) * MDP(12); (-g(1) * t32 + g(2) * t30 + t41 * t58) * MDP(18) + (g(1) * t33 - g(2) * t31 + t43 * t58) * MDP(19);];
taug = t1;
