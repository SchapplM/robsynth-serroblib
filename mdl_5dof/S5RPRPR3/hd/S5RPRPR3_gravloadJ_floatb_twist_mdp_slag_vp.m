% Calculate Gravitation load on the joints for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:32
% EndTime: 2019-12-05 17:51:33
% DurationCPUTime: 0.12s
% Computational Cost: add. (173->38), mult. (130->56), div. (0->0), fcn. (112->10), ass. (0->22)
t47 = sin(pkin(9));
t59 = g(1) * t47;
t48 = cos(pkin(9));
t49 = sin(qJ(5));
t58 = t48 * t49;
t51 = cos(qJ(5));
t57 = t48 * t51;
t46 = qJ(1) + pkin(8);
t45 = qJ(3) + t46;
t43 = sin(t45);
t44 = cos(t45);
t56 = -t43 * pkin(3) + t44 * qJ(4);
t41 = g(2) * t44 + g(3) * t43;
t54 = -t44 * pkin(3) - t43 * qJ(4);
t33 = t43 * t58 + t44 * t51;
t34 = t43 * t57 - t44 * t49;
t35 = -t43 * t51 + t44 * t58;
t36 = -t43 * t49 - t44 * t57;
t53 = (-g(2) * t35 - g(3) * t33) * MDP(18) + (-g(2) * t36 + g(3) * t34) * MDP(17) + (MDP(10) - MDP(7)) * (g(2) * t43 - g(3) * t44) + (t48 * MDP(8) - t47 * MDP(9) + MDP(6)) * t41;
t52 = cos(qJ(1));
t50 = sin(qJ(1));
t1 = [(-g(2) * t50 + g(3) * t52) * MDP(3) + (-g(2) * (-pkin(2) * cos(t46) - t52 * pkin(1) + t54) - g(3) * (-pkin(2) * sin(t46) - t50 * pkin(1) + t56)) * MDP(11) + t53 + (pkin(1) * MDP(4) + MDP(2)) * (g(2) * t52 + g(3) * t50); (-MDP(11) - MDP(4)) * g(1); (-g(2) * t54 - g(3) * t56) * MDP(11) + t53; -t41 * MDP(11); (-g(2) * t33 + g(3) * t35 + t49 * t59) * MDP(17) + (-g(2) * t34 - g(3) * t36 + t51 * t59) * MDP(18);];
taug = t1;
