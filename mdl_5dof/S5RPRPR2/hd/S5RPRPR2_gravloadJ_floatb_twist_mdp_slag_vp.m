% Calculate Gravitation load on the joints for
% S5RPRPR2
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
%   see S5RPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:38
% EndTime: 2019-12-05 17:49:38
% DurationCPUTime: 0.07s
% Computational Cost: add. (151->31), mult. (104->39), div. (0->0), fcn. (74->10), ass. (0->15)
t39 = qJ(1) + pkin(8);
t37 = qJ(3) + t39;
t33 = sin(t37);
t34 = cos(t37);
t47 = -t33 * pkin(3) + t34 * qJ(4);
t31 = g(2) * t34 + g(3) * t33;
t30 = g(2) * t33 - g(3) * t34;
t45 = -t34 * pkin(3) - t33 * qJ(4);
t38 = pkin(9) + qJ(5);
t35 = sin(t38);
t36 = cos(t38);
t44 = (MDP(10) - MDP(7)) * t30 + (t36 * MDP(17) - t35 * MDP(18) + MDP(8) * cos(pkin(9)) - MDP(9) * sin(pkin(9)) + MDP(6)) * t31;
t43 = cos(qJ(1));
t42 = sin(qJ(1));
t1 = [(-g(2) * t42 + g(3) * t43) * MDP(3) + (-g(2) * (-pkin(2) * cos(t39) - t43 * pkin(1) + t45) - g(3) * (-pkin(2) * sin(t39) - t42 * pkin(1) + t47)) * MDP(11) + t44 + (pkin(1) * MDP(4) + MDP(2)) * (g(2) * t43 + g(3) * t42); (-MDP(11) - MDP(4)) * g(1); (-g(2) * t45 - g(3) * t47) * MDP(11) + t44; -t31 * MDP(11); (-g(1) * t36 - t30 * t35) * MDP(17) + (g(1) * t35 - t30 * t36) * MDP(18);];
taug = t1;
