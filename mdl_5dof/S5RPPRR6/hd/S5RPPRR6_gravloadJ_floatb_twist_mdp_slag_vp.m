% Calculate Gravitation load on the joints for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:08
% EndTime: 2019-12-31 17:58:08
% DurationCPUTime: 0.15s
% Computational Cost: add. (117->39), mult. (118->58), div. (0->0), fcn. (102->10), ass. (0->24)
t35 = pkin(9) + qJ(4);
t31 = sin(t35);
t54 = t31 * MDP(15);
t39 = sin(qJ(5));
t41 = cos(qJ(5));
t53 = t41 * MDP(21) - t39 * MDP(22) + MDP(14);
t52 = g(3) * t31;
t36 = qJ(1) + pkin(8);
t32 = sin(t36);
t51 = t32 * t39;
t50 = t32 * t41;
t34 = cos(t36);
t49 = t34 * t39;
t48 = t34 * t41;
t45 = g(1) * t34 + g(2) * t32;
t44 = g(1) * t32 - g(2) * t34;
t42 = cos(qJ(1));
t40 = sin(qJ(1));
t33 = cos(t35);
t30 = t33 * t48 + t51;
t29 = -t33 * t49 + t50;
t28 = -t33 * t50 + t49;
t27 = t33 * t51 + t48;
t1 = [(g(1) * t42 + g(2) * t40) * MDP(3) - t45 * MDP(7) + (-g(1) * (-t40 * pkin(1) - t32 * pkin(2) + t34 * qJ(3)) - g(2) * (t42 * pkin(1) + t34 * pkin(2) + t32 * qJ(3))) * MDP(8) + (-g(1) * t28 - g(2) * t30) * MDP(21) + (-g(1) * t27 - g(2) * t29) * MDP(22) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t40 - g(2) * t42) + (t33 * MDP(14) - t54 + MDP(5) * cos(pkin(9)) - MDP(6) * sin(pkin(9))) * t44; (-MDP(4) - MDP(8)) * g(3); -t44 * MDP(8); (-t53 * t33 + t54) * g(3) + (MDP(15) * t33 + t53 * t31) * t45; (-g(1) * t29 + g(2) * t27 + t39 * t52) * MDP(21) + (g(1) * t30 - g(2) * t28 + t41 * t52) * MDP(22);];
taug = t1;
