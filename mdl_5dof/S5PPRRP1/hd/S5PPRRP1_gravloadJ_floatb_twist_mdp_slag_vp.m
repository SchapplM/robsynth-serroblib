% Calculate Gravitation load on the joints for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S5PPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:07:15
% EndTime: 2019-12-05 15:07:15
% DurationCPUTime: 0.09s
% Computational Cost: add. (80->28), mult. (106->40), div. (0->0), fcn. (89->6), ass. (0->18)
t31 = sin(pkin(7));
t32 = cos(pkin(7));
t38 = g(1) * t32 + g(2) * t31;
t30 = pkin(8) + qJ(3);
t28 = sin(t30);
t29 = cos(t30);
t24 = -g(3) * t29 + t38 * t28;
t45 = g(3) * t28;
t34 = sin(qJ(4));
t43 = t31 * t34;
t35 = cos(qJ(4));
t42 = t31 * t35;
t41 = t32 * t34;
t40 = t32 * t35;
t39 = MDP(14) + MDP(2);
t33 = -qJ(5) - pkin(6);
t27 = t35 * pkin(4) + pkin(3);
t1 = [(-MDP(1) - t39) * g(3); t39 * (-g(1) * t31 + g(2) * t32); (-g(3) * (t29 * t27 - t28 * t33) + t38 * (t27 * t28 + t29 * t33)) * MDP(14) + (MDP(5) - MDP(13)) * (t38 * t29 + t45) + (MDP(11) * t35 - MDP(12) * t34 + MDP(4)) * t24; (-g(1) * (-t29 * t40 - t43) - g(2) * (-t29 * t42 + t41) + t35 * t45) * MDP(12) + (pkin(4) * MDP(14) + MDP(11)) * (-g(1) * (-t29 * t41 + t42) - g(2) * (-t29 * t43 - t40) + t34 * t45); -t24 * MDP(14);];
taug = t1;
