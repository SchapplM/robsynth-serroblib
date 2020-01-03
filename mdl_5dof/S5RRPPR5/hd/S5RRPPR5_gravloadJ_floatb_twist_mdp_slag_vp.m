% Calculate Gravitation load on the joints for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:29:58
% EndTime: 2019-12-31 19:29:59
% DurationCPUTime: 0.20s
% Computational Cost: add. (149->50), mult. (200->71), div. (0->0), fcn. (181->8), ass. (0->25)
t52 = qJ(2) + pkin(8);
t48 = sin(t52);
t49 = cos(t52);
t54 = sin(qJ(5));
t57 = cos(qJ(5));
t71 = -t48 * t57 + t49 * t54;
t56 = sin(qJ(1));
t59 = cos(qJ(1));
t45 = g(1) * t59 + g(2) * t56;
t37 = t71 * t56;
t64 = t48 * t54 + t49 * t57;
t38 = t64 * t56;
t39 = t71 * t59;
t40 = t64 * t59;
t70 = (g(1) * t39 + g(2) * t37 + g(3) * t64) * MDP(22) + (g(1) * t40 + g(2) * t38 - g(3) * t71) * MDP(23);
t44 = g(1) * t56 - g(2) * t59;
t65 = t49 * pkin(3) + t48 * qJ(4);
t58 = cos(qJ(2));
t55 = sin(qJ(2));
t53 = -qJ(3) - pkin(6);
t50 = t58 * pkin(2);
t47 = t50 + pkin(1);
t46 = t59 * t47;
t36 = -g(3) * t49 + t45 * t48;
t1 = [(-g(1) * (-t56 * t47 - t59 * t53) - g(2) * (-t56 * t53 + t46)) * MDP(12) + (-g(2) * t46 + (g(1) * t53 - g(2) * t65) * t59 + (-g(1) * (-t47 - t65) + g(2) * t53) * t56) * MDP(16) + (g(1) * t38 - g(2) * t40) * MDP(22) + (-g(1) * t37 + g(2) * t39) * MDP(23) + (MDP(3) - MDP(11) - MDP(14)) * t45 + (-t55 * MDP(10) + t49 * MDP(13) + t48 * MDP(15) + t58 * MDP(9) + MDP(2)) * t44; (g(3) * t55 + t45 * t58) * MDP(10) + t36 * MDP(13) + (-g(3) * t48 - t45 * t49) * MDP(15) + (-g(3) * (t50 + t65) + t45 * (pkin(2) * t55 + pkin(3) * t48 - qJ(4) * t49)) * MDP(16) - t70 + (pkin(2) * MDP(12) + MDP(9)) * (-g(3) * t58 + t45 * t55); (-MDP(12) - MDP(16)) * t44; -t36 * MDP(16); t70;];
taug = t1;
