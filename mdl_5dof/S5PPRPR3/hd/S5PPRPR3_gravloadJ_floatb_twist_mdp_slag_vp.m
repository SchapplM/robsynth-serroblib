% Calculate Gravitation load on the joints for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S5PPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:23
% EndTime: 2019-12-05 15:05:23
% DurationCPUTime: 0.14s
% Computational Cost: add. (66->35), mult. (119->64), div. (0->0), fcn. (121->10), ass. (0->25)
t36 = sin(pkin(8));
t56 = g(3) * t36;
t40 = sin(qJ(5));
t55 = t36 * t40;
t42 = cos(qJ(5));
t54 = t36 * t42;
t37 = sin(pkin(7));
t38 = cos(pkin(8));
t53 = t37 * t38;
t41 = sin(qJ(3));
t52 = t37 * t41;
t43 = cos(qJ(3));
t51 = t37 * t43;
t35 = qJ(3) + pkin(9);
t33 = sin(t35);
t39 = cos(pkin(7));
t50 = t39 * t33;
t34 = cos(t35);
t49 = t39 * t34;
t48 = t39 * t41;
t47 = t39 * t43;
t46 = MDP(2) + MDP(6);
t31 = t37 * t33 + t38 * t49;
t29 = t34 * t53 - t50;
t1 = [(-MDP(1) - t46) * g(3); t46 * (-g(1) * t37 + g(2) * t39); (-g(1) * (-t38 * t47 - t52) - g(2) * (-t38 * t51 + t48) + t43 * t56) * MDP(5) + (MDP(12) * t42 - MDP(13) * t40) * (-g(1) * (t37 * t34 - t38 * t50) - g(2) * (-t33 * t53 - t49) + t33 * t56) + (pkin(3) * MDP(6) + MDP(4)) * (-g(1) * (-t38 * t48 + t51) - g(2) * (-t38 * t52 - t47) + t41 * t56); (g(3) * t38 + (-g(1) * t39 - g(2) * t37) * t36) * MDP(6); (-g(1) * (-t31 * t40 + t39 * t54) - g(2) * (-t29 * t40 + t37 * t54) - g(3) * (-t34 * t55 - t38 * t42)) * MDP(12) + (-g(1) * (-t31 * t42 - t39 * t55) - g(2) * (-t29 * t42 - t37 * t55) - g(3) * (-t34 * t54 + t38 * t40)) * MDP(13);];
taug = t1;
