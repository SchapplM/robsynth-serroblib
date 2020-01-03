% Calculate Gravitation load on the joints for
% S5RPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RPRRR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:08:00
% EndTime: 2019-12-31 19:08:00
% DurationCPUTime: 0.14s
% Computational Cost: add. (153->39), mult. (161->56), div. (0->0), fcn. (142->10), ass. (0->23)
t44 = pkin(9) + qJ(3);
t43 = qJ(4) + t44;
t39 = sin(t43);
t59 = g(3) * t39;
t47 = sin(qJ(5));
t48 = sin(qJ(1));
t57 = t48 * t47;
t49 = cos(qJ(5));
t56 = t48 * t49;
t50 = cos(qJ(1));
t55 = t50 * t47;
t54 = t50 * t49;
t40 = cos(t43);
t52 = g(1) * t50 + g(2) * t48;
t53 = (t40 * t52 + t59) * MDP(21) + (MDP(27) * t49 - t47 * MDP(28) + MDP(20)) * (-g(3) * t40 + t39 * t52);
t37 = g(1) * t48 - g(2) * t50;
t42 = cos(t44);
t41 = sin(t44);
t36 = t40 * t54 + t57;
t35 = -t40 * t55 + t56;
t34 = -t40 * t56 + t55;
t33 = t40 * t57 + t54;
t1 = [(-g(1) * (-t48 * pkin(1) + qJ(2) * t50) - g(2) * (pkin(1) * t50 + t48 * qJ(2))) * MDP(7) + (-g(1) * t34 - g(2) * t36) * MDP(27) + (-g(1) * t33 - g(2) * t35) * MDP(28) + (MDP(3) - MDP(6)) * t52 + (t42 * MDP(13) - t41 * MDP(14) + MDP(20) * t40 - MDP(21) * t39 + MDP(4) * cos(pkin(9)) - MDP(5) * sin(pkin(9)) + MDP(2)) * t37; -t37 * MDP(7); (-g(3) * t42 + t41 * t52) * MDP(13) + (g(3) * t41 + t42 * t52) * MDP(14) + t53; t53; (-g(1) * t35 + g(2) * t33 + t47 * t59) * MDP(27) + (g(1) * t36 - g(2) * t34 + t49 * t59) * MDP(28);];
taug = t1;
