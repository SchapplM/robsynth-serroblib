% Calculate Gravitation load on the joints for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:42
% EndTime: 2019-12-05 15:24:43
% DurationCPUTime: 0.21s
% Computational Cost: add. (94->33), mult. (119->49), div. (0->0), fcn. (100->10), ass. (0->19)
t34 = sin(pkin(7));
t36 = cos(pkin(7));
t43 = g(1) * t36 + g(2) * t34;
t32 = qJ(2) + pkin(8);
t28 = sin(t32);
t30 = cos(t32);
t41 = -g(3) * t30 + t43 * t28;
t50 = g(3) * t28;
t31 = pkin(9) + qJ(5);
t27 = sin(t31);
t48 = t34 * t27;
t29 = cos(t31);
t47 = t34 * t29;
t46 = t36 * t27;
t45 = t36 * t29;
t44 = MDP(5) + MDP(9);
t38 = cos(qJ(2));
t37 = sin(qJ(2));
t1 = [(-MDP(1) - t44) * g(3); (g(3) * t37 + t43 * t38) * MDP(4) + (-t43 * t30 - t50) * MDP(8) + (-g(3) * (t38 * pkin(2) + t30 * pkin(3) + t28 * qJ(4)) + t43 * (pkin(2) * t37 + pkin(3) * t28 - qJ(4) * t30)) * MDP(9) + (pkin(2) * MDP(5) + MDP(3)) * (-g(3) * t38 + t43 * t37) + (MDP(15) * t29 - MDP(16) * t27 + MDP(6) * cos(pkin(9)) - MDP(7) * sin(pkin(9))) * t41; t44 * (-g(1) * t34 + g(2) * t36); -t41 * MDP(9); (-g(1) * (-t30 * t46 + t47) - g(2) * (-t30 * t48 - t45) + t27 * t50) * MDP(15) + (-g(1) * (-t30 * t45 - t48) - g(2) * (-t30 * t47 + t46) + t29 * t50) * MDP(16);];
taug = t1;
