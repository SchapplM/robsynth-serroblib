% Calculate Gravitation load on the joints for
% S5PPRRR2
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
%   see S5PPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PPRRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:49
% EndTime: 2019-12-05 15:14:49
% DurationCPUTime: 0.13s
% Computational Cost: add. (109->31), mult. (119->49), div. (0->0), fcn. (116->8), ass. (0->22)
t33 = qJ(4) + qJ(5);
t30 = sin(t33);
t31 = cos(t33);
t36 = sin(qJ(4));
t37 = cos(qJ(4));
t53 = t37 * MDP(11) - t36 * MDP(12) + t31 * MDP(18) - t30 * MDP(19) + MDP(4);
t32 = pkin(9) + qJ(3);
t28 = sin(t32);
t52 = g(3) * t28;
t34 = sin(pkin(8));
t51 = t34 * t30;
t50 = t34 * t31;
t49 = t34 * t36;
t48 = t34 * t37;
t35 = cos(pkin(8));
t47 = t35 * t30;
t46 = t35 * t31;
t45 = t35 * t36;
t44 = t35 * t37;
t29 = cos(t32);
t43 = (-g(1) * (-t29 * t47 + t50) - g(2) * (-t29 * t51 - t46) + t30 * t52) * MDP(18) + (-g(1) * (-t29 * t46 - t51) - g(2) * (-t29 * t50 + t47) + t31 * t52) * MDP(19);
t1 = [(-MDP(1) - MDP(2)) * g(3); (-g(1) * t34 + g(2) * t35) * MDP(2); (t28 * MDP(5) - t53 * t29) * g(3) + (MDP(5) * t29 + t53 * t28) * (g(1) * t35 + g(2) * t34); (-g(1) * (-t29 * t45 + t48) - g(2) * (-t29 * t49 - t44) + t36 * t52) * MDP(11) + (-g(1) * (-t29 * t44 - t49) - g(2) * (-t29 * t48 + t45) + t37 * t52) * MDP(12) + t43; t43;];
taug = t1;
