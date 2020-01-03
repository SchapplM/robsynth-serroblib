% Calculate Gravitation load on the joints for
% S4RRPR9
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
%   see S4RRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:05
% EndTime: 2019-12-31 17:10:05
% DurationCPUTime: 0.20s
% Computational Cost: add. (94->45), mult. (168->70), div. (0->0), fcn. (154->8), ass. (0->25)
t44 = sin(qJ(1));
t46 = cos(qJ(1));
t52 = g(1) * t46 + g(2) * t44;
t62 = MDP(10) - MDP(13);
t43 = sin(qJ(2));
t45 = cos(qJ(2));
t35 = -g(3) * t45 + t52 * t43;
t59 = g(3) * t43;
t57 = t44 * t45;
t40 = pkin(7) + qJ(4);
t38 = sin(t40);
t56 = t46 * t38;
t39 = cos(t40);
t55 = t46 * t39;
t41 = sin(pkin(7));
t54 = t46 * t41;
t42 = cos(pkin(7));
t53 = t46 * t42;
t50 = t45 * pkin(2) + t43 * qJ(3);
t48 = pkin(1) + t50;
t34 = t44 * t38 + t45 * t55;
t33 = t44 * t39 - t45 * t56;
t32 = -t39 * t57 + t56;
t31 = t38 * t57 + t55;
t1 = [t52 * MDP(3) + (-g(1) * (-t42 * t57 + t54) - g(2) * (t44 * t41 + t45 * t53)) * MDP(11) + (-g(1) * (t41 * t57 + t53) - g(2) * (t44 * t42 - t45 * t54)) * MDP(12) + ((-g(1) * pkin(5) - g(2) * t48) * t46 + (-g(2) * pkin(5) + g(1) * t48) * t44) * MDP(14) + (-g(1) * t32 - g(2) * t34) * MDP(20) + (-g(1) * t31 - g(2) * t33) * MDP(21) + (t45 * MDP(9) - t62 * t43 + MDP(2)) * (g(1) * t44 - g(2) * t46); (-g(3) * t50 + t52 * (pkin(2) * t43 - qJ(3) * t45)) * MDP(14) + t62 * (t52 * t45 + t59) + (MDP(11) * t42 - MDP(12) * t41 + MDP(20) * t39 - MDP(21) * t38 + MDP(9)) * t35; -t35 * MDP(14); (-g(1) * t33 + g(2) * t31 + t38 * t59) * MDP(20) + (g(1) * t34 - g(2) * t32 + t39 * t59) * MDP(21);];
taug = t1;
