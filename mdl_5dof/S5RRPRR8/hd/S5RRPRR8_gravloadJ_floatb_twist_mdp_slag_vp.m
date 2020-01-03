% Calculate Gravitation load on the joints for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:13
% EndTime: 2019-12-31 20:18:14
% DurationCPUTime: 0.12s
% Computational Cost: add. (146->39), mult. (162->56), div. (0->0), fcn. (141->8), ass. (0->24)
t42 = qJ(2) + pkin(9) + qJ(4);
t39 = sin(t42);
t59 = g(3) * t39;
t44 = sin(qJ(5));
t46 = sin(qJ(1));
t57 = t46 * t44;
t47 = cos(qJ(5));
t56 = t46 * t47;
t49 = cos(qJ(1));
t55 = t49 * t44;
t54 = t49 * t47;
t40 = cos(t42);
t52 = g(1) * t49 + g(2) * t46;
t53 = (t52 * t40 + t59) * MDP(19) + (t47 * MDP(25) - t44 * MDP(26) + MDP(18)) * (-g(3) * t40 + t52 * t39);
t37 = g(1) * t46 - g(2) * t49;
t48 = cos(qJ(2));
t45 = sin(qJ(2));
t43 = -qJ(3) - pkin(6);
t41 = t48 * pkin(2) + pkin(1);
t36 = t40 * t54 + t57;
t35 = -t40 * t55 + t56;
t34 = -t40 * t56 + t55;
t33 = t40 * t57 + t54;
t1 = [(-g(1) * (-t46 * t41 - t49 * t43) - g(2) * (t49 * t41 - t46 * t43)) * MDP(12) + (-g(1) * t34 - g(2) * t36) * MDP(25) + (-g(1) * t33 - g(2) * t35) * MDP(26) + (MDP(3) - MDP(11)) * t52 + (-t45 * MDP(10) + MDP(18) * t40 - MDP(19) * t39 + t48 * MDP(9) + MDP(2)) * t37; (g(3) * t45 + t52 * t48) * MDP(10) + t53 + (MDP(12) * pkin(2) + MDP(9)) * (-g(3) * t48 + t52 * t45); -t37 * MDP(12); t53; (-g(1) * t35 + g(2) * t33 + t44 * t59) * MDP(25) + (g(1) * t36 - g(2) * t34 + t47 * t59) * MDP(26);];
taug = t1;
