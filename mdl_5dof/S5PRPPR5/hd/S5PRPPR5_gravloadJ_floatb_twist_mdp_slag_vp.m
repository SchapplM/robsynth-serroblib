% Calculate Gravitation load on the joints for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:38
% EndTime: 2019-12-31 17:38:38
% DurationCPUTime: 0.20s
% Computational Cost: add. (73->39), mult. (175->62), div. (0->0), fcn. (173->8), ass. (0->23)
t47 = sin(pkin(7));
t49 = cos(pkin(7));
t56 = g(1) * t49 + g(2) * t47;
t51 = sin(qJ(2));
t65 = t56 * t51;
t63 = pkin(2) * t51;
t46 = sin(pkin(8));
t48 = cos(pkin(8));
t53 = cos(qJ(2));
t55 = t53 * t46 - t51 * t48;
t60 = g(3) * t55;
t59 = t53 * pkin(2) + t51 * qJ(3);
t58 = qJ(3) * t53;
t57 = -MDP(10) - MDP(7);
t40 = t51 * t46 + t53 * t48;
t52 = cos(qJ(5));
t50 = sin(qJ(5));
t42 = t49 * t58;
t41 = t47 * t58;
t38 = t40 * t49;
t36 = t40 * t47;
t33 = -g(3) * t53 + t65;
t1 = [(-MDP(1) + t57) * g(3); (-g(1) * (-t49 * t63 + t42) - g(2) * (-t47 * t63 + t41) - g(3) * t59) * MDP(7) + (-g(1) * t38 - g(2) * t36 + t60) * MDP(9) + (-g(1) * t42 - g(2) * t41 - g(3) * (t53 * pkin(3) + t59) + (pkin(2) + pkin(3)) * t65) * MDP(10) + (MDP(4) - MDP(6)) * (g(3) * t51 + t56 * t53) + (MDP(3) + MDP(5)) * t33 + (-MDP(16) * t52 + MDP(17) * t50 - MDP(8)) * (g(3) * t40 + t56 * t55); t57 * t33; (g(1) * t47 - g(2) * t49) * MDP(10); (-g(1) * (-t38 * t50 - t47 * t52) - g(2) * (-t36 * t50 + t49 * t52) - t50 * t60) * MDP(16) + (-g(1) * (-t38 * t52 + t47 * t50) - g(2) * (-t36 * t52 - t49 * t50) - t52 * t60) * MDP(17);];
taug = t1;
