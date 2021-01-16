% Calculate Gravitation load on the joints for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:56:39
% EndTime: 2021-01-15 10:56:40
% DurationCPUTime: 0.14s
% Computational Cost: add. (85->38), mult. (132->54), div. (0->0), fcn. (115->8), ass. (0->25)
t37 = sin(qJ(4));
t40 = cos(qJ(4));
t53 = t40 * MDP(20) - t37 * MDP(21) + MDP(11);
t35 = qJ(2) + pkin(7);
t33 = sin(t35);
t38 = sin(qJ(2));
t52 = t38 * MDP(10) + t33 * MDP(12);
t51 = g(3) * t33;
t39 = sin(qJ(1));
t50 = t39 * t37;
t49 = t39 * t40;
t42 = cos(qJ(1));
t48 = t42 * t37;
t47 = t42 * t40;
t44 = g(1) * t42 + g(2) * t39;
t30 = g(1) * t39 - g(2) * t42;
t41 = cos(qJ(2));
t36 = -qJ(3) - pkin(5);
t34 = cos(t35);
t32 = t41 * pkin(2) + pkin(1);
t29 = t34 * t47 + t50;
t28 = -t34 * t48 + t49;
t27 = -t34 * t49 + t48;
t26 = t34 * t50 + t47;
t1 = [(-g(1) * (-t39 * t32 - t36 * t42) - g(2) * (t42 * t32 - t39 * t36)) * MDP(14) + (-g(1) * t27 - g(2) * t29) * MDP(20) + (-g(1) * t26 - g(2) * t28) * MDP(21) + (MDP(3) - MDP(13)) * t44 + (t34 * MDP(11) + t41 * MDP(9) + MDP(2) - t52) * t30; (-t53 * t34 + t52) * g(3) + (t41 * MDP(10) + MDP(12) * t34 + t53 * t33) * t44 + (MDP(14) * pkin(2) + MDP(9)) * (-g(3) * t41 + t44 * t38); -t30 * MDP(14); (-g(1) * t28 + g(2) * t26 + t37 * t51) * MDP(20) + (g(1) * t29 - g(2) * t27 + t40 * t51) * MDP(21);];
taug = t1;
