% Calculate Gravitation load on the joints for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:29:06
% EndTime: 2019-12-05 17:29:07
% DurationCPUTime: 0.21s
% Computational Cost: add. (129->50), mult. (128->79), div. (0->0), fcn. (113->10), ass. (0->29)
t47 = sin(pkin(8));
t64 = g(1) * t47;
t45 = qJ(1) + pkin(7);
t43 = cos(t45);
t63 = g(2) * t43;
t51 = cos(qJ(1));
t62 = t51 * pkin(1);
t41 = sin(t45);
t49 = cos(pkin(8));
t61 = t41 * t49;
t60 = t43 * t49;
t46 = sin(pkin(9));
t59 = t46 * t49;
t48 = cos(pkin(9));
t58 = t48 * t49;
t57 = MDP(12) + MDP(8);
t50 = sin(qJ(1));
t56 = -t50 * pkin(1) + t43 * qJ(3);
t55 = g(3) * t41 + t63;
t54 = g(2) * t41 - g(3) * t43;
t52 = pkin(3) * t49 + qJ(4) * t47 + pkin(2);
t44 = pkin(9) + qJ(5);
t42 = cos(t44);
t40 = sin(t44);
t36 = -t41 * t40 - t42 * t60;
t35 = t40 * t60 - t41 * t42;
t34 = -t43 * t40 + t42 * t61;
t33 = t40 * t61 + t43 * t42;
t1 = [(-g(2) * t50 + g(3) * t51) * MDP(3) + t54 * MDP(7) + (-g(2) * (-t43 * pkin(2) - t41 * qJ(3) - t62) - g(3) * (-t41 * pkin(2) + t56)) * MDP(8) + (-g(2) * (-t41 * t46 - t43 * t58) - g(3) * (-t41 * t58 + t43 * t46)) * MDP(9) + (-g(2) * (-t41 * t48 + t43 * t59) - g(3) * (t41 * t59 + t43 * t48)) * MDP(10) + (g(2) * t62 - g(3) * t56 + t52 * t63 + (g(2) * qJ(3) + g(3) * t52) * t41) * MDP(12) + (-g(2) * t36 + g(3) * t34) * MDP(18) + (-g(2) * t35 - g(3) * t33) * MDP(19) + (pkin(1) * MDP(4) + MDP(2)) * (g(2) * t51 + g(3) * t50) + (t49 * MDP(5) + (-MDP(6) + MDP(11)) * t47) * t55; (-MDP(4) - t57) * g(1); -t57 * t55; (g(1) * t49 + t54 * t47) * MDP(12); (-g(2) * t33 + g(3) * t35 + t40 * t64) * MDP(18) + (-g(2) * t34 - g(3) * t36 + t42 * t64) * MDP(19);];
taug = t1;
