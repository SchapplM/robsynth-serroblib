% Calculate Gravitation load on the joints for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPPPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:45
% EndTime: 2022-01-20 09:12:46
% DurationCPUTime: 0.16s
% Computational Cost: add. (121->47), mult. (118->76), div. (0->0), fcn. (105->10), ass. (0->29)
t48 = qJ(1) + pkin(7);
t43 = sin(t48);
t67 = g(1) * t43;
t50 = sin(pkin(8));
t66 = g(3) * t50;
t52 = cos(pkin(8));
t65 = t43 * t52;
t45 = cos(t48);
t64 = t45 * t52;
t49 = sin(pkin(9));
t63 = t49 * t52;
t51 = cos(pkin(9));
t62 = t51 * t52;
t61 = MDP(10) + MDP(7);
t54 = cos(qJ(1));
t60 = t54 * pkin(1) + t45 * pkin(2) + t43 * qJ(3);
t53 = sin(qJ(1));
t59 = -t53 * pkin(1) + t45 * qJ(3);
t58 = -g(1) * t45 - g(2) * t43;
t57 = -g(2) * t45 + t67;
t55 = pkin(3) * t52 + qJ(4) * t50;
t47 = pkin(9) + qJ(5);
t44 = cos(t47);
t42 = sin(t47);
t37 = t43 * t42 + t44 * t64;
t36 = -t42 * t64 + t43 * t44;
t35 = t45 * t42 - t44 * t65;
t34 = t42 * t65 + t45 * t44;
t1 = [(g(1) * t54 + g(2) * t53) * MDP(3) + t57 * t52 * MDP(5) + t58 * MDP(6) + (-g(1) * (-t43 * pkin(2) + t59) - g(2) * t60) * MDP(7) + (-g(1) * (-t43 * t62 + t45 * t49) - g(2) * (t43 * t49 + t45 * t62)) * MDP(8) + (-g(1) * (t43 * t63 + t45 * t51) - g(2) * (t43 * t51 - t45 * t63)) * MDP(9) + (-g(1) * t59 - g(2) * (t55 * t45 + t60) - (-pkin(2) - t55) * t67) * MDP(10) + (-g(1) * t35 - g(2) * t37) * MDP(16) + (-g(1) * t34 - g(2) * t36) * MDP(17) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t53 - g(2) * t54); (-MDP(4) - t61) * g(3); -t61 * t57; (g(3) * t52 + t58 * t50) * MDP(10); (-g(1) * t36 + g(2) * t34 + t42 * t66) * MDP(16) + (g(1) * t37 - g(2) * t35 + t44 * t66) * MDP(17);];
taug = t1;
