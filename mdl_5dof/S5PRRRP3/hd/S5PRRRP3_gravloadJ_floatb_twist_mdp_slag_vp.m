% Calculate Gravitation load on the joints for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:44:13
% EndTime: 2019-12-05 16:44:13
% DurationCPUTime: 0.07s
% Computational Cost: add. (112->30), mult. (102->40), div. (0->0), fcn. (75->6), ass. (0->16)
t33 = qJ(3) + qJ(4);
t28 = sin(t33);
t29 = cos(t33);
t32 = pkin(8) + qJ(2);
t26 = sin(t32);
t27 = cos(t32);
t37 = g(1) * t27 + g(2) * t26;
t36 = -g(3) * t29 + t37 * t28;
t39 = t36 * MDP(17) + (g(3) * t28 + t37 * t29) * MDP(18);
t35 = cos(qJ(3));
t38 = t35 * pkin(3) + pkin(4) * t29;
t21 = g(1) * t26 - g(2) * t27;
t34 = sin(qJ(3));
t31 = -qJ(5) - pkin(7) - pkin(6);
t23 = pkin(2) + t38;
t1 = [(-MDP(1) - MDP(20)) * g(3); (-g(1) * (-t26 * t23 - t27 * t31) - g(2) * (t27 * t23 - t26 * t31)) * MDP(20) + (MDP(4) - MDP(19)) * t37 + (t35 * MDP(10) - t34 * MDP(11) + MDP(17) * t29 - MDP(18) * t28 + MDP(3)) * t21; (-g(3) * t35 + t37 * t34) * MDP(10) + (g(3) * t34 + t37 * t35) * MDP(11) + (-g(3) * t38 - t37 * (-t34 * pkin(3) - pkin(4) * t28)) * MDP(20) + t39; t36 * MDP(20) * pkin(4) + t39; -t21 * MDP(20);];
taug = t1;
