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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:36:21
% EndTime: 2020-01-03 11:36:22
% DurationCPUTime: 0.09s
% Computational Cost: add. (173->37), mult. (130->56), div. (0->0), fcn. (112->10), ass. (0->22)
t51 = sin(pkin(9));
t63 = g(1) * t51;
t52 = cos(pkin(9));
t53 = sin(qJ(5));
t62 = t52 * t53;
t55 = cos(qJ(5));
t61 = t52 * t55;
t50 = qJ(1) + pkin(8);
t49 = qJ(3) + t50;
t47 = sin(t49);
t48 = cos(t49);
t60 = t48 * pkin(3) + t47 * qJ(4);
t59 = t47 * pkin(3) - t48 * qJ(4);
t43 = g(2) * t48 + g(3) * t47;
t35 = -t47 * t62 - t48 * t55;
t36 = t47 * t61 - t48 * t53;
t37 = -t47 * t55 + t48 * t62;
t38 = t47 * t53 + t48 * t61;
t57 = (g(2) * t37 - g(3) * t35) * MDP(18) + (-g(2) * t38 - g(3) * t36) * MDP(17) + (MDP(7) - MDP(10)) * (g(2) * t47 - g(3) * t48) + (-t52 * MDP(8) + t51 * MDP(9) - MDP(6)) * t43;
t56 = cos(qJ(1));
t54 = sin(qJ(1));
t1 = [(g(2) * t54 - g(3) * t56) * MDP(3) + (-g(2) * (pkin(2) * cos(t50) + t56 * pkin(1) + t60) - g(3) * (pkin(2) * sin(t50) + t54 * pkin(1) + t59)) * MDP(11) + t57 + (pkin(1) * MDP(4) + MDP(2)) * (-g(2) * t56 - g(3) * t54); (-MDP(11) - MDP(4)) * g(1); (-g(2) * t60 - g(3) * t59) * MDP(11) + t57; t43 * MDP(11); (-g(2) * t35 - g(3) * t37 + t53 * t63) * MDP(17) + (g(2) * t36 - g(3) * t38 + t55 * t63) * MDP(18);];
taug = t1;
