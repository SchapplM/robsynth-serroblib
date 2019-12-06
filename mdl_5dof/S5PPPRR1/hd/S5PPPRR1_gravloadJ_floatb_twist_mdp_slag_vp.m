% Calculate Gravitation load on the joints for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPPRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S5PPPRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:58:07
% EndTime: 2019-12-05 14:58:08
% DurationCPUTime: 0.07s
% Computational Cost: add. (71->28), mult. (103->51), div. (0->0), fcn. (109->8), ass. (0->19)
t29 = sin(pkin(8));
t42 = g(3) * t29;
t33 = sin(qJ(5));
t41 = t29 * t33;
t34 = cos(qJ(5));
t40 = t29 * t34;
t30 = sin(pkin(7));
t31 = cos(pkin(8));
t39 = t30 * t31;
t28 = pkin(9) + qJ(4);
t26 = sin(t28);
t32 = cos(pkin(7));
t38 = t32 * t26;
t27 = cos(t28);
t37 = t32 * t27;
t36 = MDP(2) + MDP(3);
t24 = t30 * t26 + t31 * t37;
t22 = t27 * t39 - t38;
t1 = [(-MDP(1) - t36) * g(3); t36 * (-g(1) * t30 + g(2) * t32); (g(3) * t31 + (-g(1) * t32 - g(2) * t30) * t29) * MDP(3); (g(1) * t24 + g(2) * t22 + t27 * t42) * MDP(6) + (MDP(12) * t34 - MDP(13) * t33 + MDP(5)) * (-g(1) * (t30 * t27 - t31 * t38) - g(2) * (-t26 * t39 - t37) + t26 * t42); (-g(1) * (-t24 * t33 + t32 * t40) - g(2) * (-t22 * t33 + t30 * t40) - g(3) * (-t27 * t41 - t31 * t34)) * MDP(12) + (-g(1) * (-t24 * t34 - t32 * t41) - g(2) * (-t22 * t34 - t30 * t41) - g(3) * (-t27 * t40 + t31 * t33)) * MDP(13);];
taug = t1;
