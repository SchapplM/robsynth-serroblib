% Calculate Gravitation load on the joints for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:04
% EndTime: 2019-12-05 15:31:05
% DurationCPUTime: 0.17s
% Computational Cost: add. (110->39), mult. (123->58), div. (0->0), fcn. (105->6), ass. (0->22)
t48 = sin(pkin(8));
t49 = cos(pkin(8));
t52 = cos(qJ(4));
t65 = -t48 * (-qJ(5) - pkin(6)) + (t52 * pkin(4) + pkin(3)) * t49;
t61 = g(3) * t48;
t47 = pkin(7) + qJ(2);
t46 = cos(t47);
t51 = sin(qJ(4));
t59 = t46 * t51;
t57 = t49 * t51;
t56 = t49 * t52;
t45 = sin(t47);
t55 = t46 * pkin(2) + t45 * qJ(3);
t54 = -MDP(17) - MDP(8);
t39 = g(1) * t46 + g(2) * t45;
t38 = g(1) * t45 - g(2) * t46;
t35 = t45 * t52 - t46 * t57;
t33 = t45 * t57 + t46 * t52;
t41 = t46 * qJ(3);
t36 = t45 * t51 + t46 * t56;
t34 = -t45 * t56 + t59;
t1 = [(-MDP(1) + t54) * g(3); (-g(1) * (-t45 * pkin(2) + t41) - g(2) * t55) * MDP(8) + (-g(1) * t34 - g(2) * t36) * MDP(14) + (-g(1) * t33 - g(2) * t35) * MDP(15) + (-g(1) * (pkin(4) * t59 + t41) - g(2) * (t65 * t46 + t55) + (-g(1) * (-pkin(2) - t65) - g(2) * pkin(4) * t51) * t45) * MDP(17) + (MDP(4) - MDP(7)) * t39 + (t49 * MDP(5) + MDP(3) + (-MDP(6) + MDP(16)) * t48) * t38; t54 * t38; (g(1) * t36 - g(2) * t34 + t52 * t61) * MDP(15) + (pkin(4) * MDP(17) + MDP(14)) * (-g(1) * t35 + g(2) * t33 + t51 * t61); (g(3) * t49 - t39 * t48) * MDP(17);];
taug = t1;
