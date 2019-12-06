% Calculate Gravitation load on the joints for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:35
% EndTime: 2019-12-05 16:17:35
% DurationCPUTime: 0.07s
% Computational Cost: add. (173->35), mult. (123->53), div. (0->0), fcn. (108->8), ass. (0->22)
t51 = sin(pkin(9));
t60 = g(3) * t51;
t52 = cos(pkin(9));
t53 = sin(qJ(5));
t59 = t52 * t53;
t54 = cos(qJ(5));
t58 = t52 * t54;
t50 = pkin(8) + qJ(2);
t49 = qJ(3) + t50;
t45 = sin(t49);
t46 = cos(t49);
t57 = t46 * pkin(3) + t45 * qJ(4);
t56 = -t45 * pkin(3) + t46 * qJ(4);
t40 = g(1) * t45 - g(2) * t46;
t33 = t45 * t59 + t46 * t54;
t34 = -t45 * t58 + t46 * t53;
t35 = t45 * t54 - t46 * t59;
t36 = t45 * t53 + t46 * t58;
t55 = (-g(1) * t33 - g(2) * t35) * MDP(18) + (-g(1) * t34 - g(2) * t36) * MDP(17) + (MDP(7) - MDP(10)) * (g(1) * t46 + g(2) * t45) + (t52 * MDP(8) - t51 * MDP(9) + MDP(6)) * t40;
t48 = cos(t50);
t47 = sin(t50);
t1 = [(-MDP(1) - MDP(11)) * g(3); (g(1) * t47 - g(2) * t48) * MDP(3) + (g(1) * t48 + g(2) * t47) * MDP(4) + (-g(1) * (-pkin(2) * t47 + t56) - g(2) * (pkin(2) * t48 + t57)) * MDP(11) + t55; (-g(1) * t56 - g(2) * t57) * MDP(11) + t55; -t40 * MDP(11); (-g(1) * t35 + g(2) * t33 + t53 * t60) * MDP(17) + (g(1) * t36 - g(2) * t34 + t54 * t60) * MDP(18);];
taug = t1;
