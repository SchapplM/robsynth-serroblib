% Calculate Gravitation load on the joints for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4PRRR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:48
% EndTime: 2019-12-31 16:36:49
% DurationCPUTime: 0.18s
% Computational Cost: add. (96->48), mult. (251->91), div. (0->0), fcn. (300->10), ass. (0->28)
t46 = sin(pkin(4));
t66 = g(3) * t46;
t50 = sin(qJ(3));
t65 = t46 * t50;
t51 = sin(qJ(2));
t64 = t46 * t51;
t53 = cos(qJ(3));
t63 = t46 * t53;
t54 = cos(qJ(2));
t62 = t46 * t54;
t48 = cos(pkin(4));
t61 = t48 * t51;
t60 = t48 * t54;
t49 = sin(qJ(4));
t59 = t49 * t53;
t52 = cos(qJ(4));
t58 = t52 * t53;
t57 = t53 * t54;
t47 = cos(pkin(8));
t45 = sin(pkin(8));
t43 = t48 * t50 + t51 * t63;
t41 = -t45 * t61 + t47 * t54;
t40 = t45 * t60 + t47 * t51;
t39 = t45 * t54 + t47 * t61;
t38 = t45 * t51 - t47 * t60;
t37 = t41 * t53 + t45 * t65;
t35 = t39 * t53 - t47 * t65;
t1 = [-g(3) * MDP(1); (g(1) * t41 + g(2) * t39 + g(3) * t64) * MDP(4) + (-g(1) * (-t40 * t58 + t41 * t49) - g(2) * (-t38 * t58 + t39 * t49) - (t49 * t51 + t52 * t57) * t66) * MDP(17) + (-g(1) * (t40 * t59 + t41 * t52) - g(2) * (t38 * t59 + t39 * t52) - (-t49 * t57 + t51 * t52) * t66) * MDP(18) + (-t53 * MDP(10) + MDP(11) * t50 - MDP(3)) * (-g(1) * t40 - g(2) * t38 + g(3) * t62); (g(1) * t37 + g(2) * t35 + g(3) * t43) * MDP(11) + (-MDP(17) * t52 + MDP(18) * t49 - MDP(10)) * (g(1) * (-t41 * t50 + t45 * t63) + g(2) * (-t39 * t50 - t47 * t63) + g(3) * (t48 * t53 - t50 * t64)); (-g(1) * (-t37 * t49 + t40 * t52) - g(2) * (-t35 * t49 + t38 * t52) - g(3) * (-t43 * t49 - t52 * t62)) * MDP(17) + (-g(1) * (-t37 * t52 - t40 * t49) - g(2) * (-t35 * t52 - t38 * t49) - g(3) * (-t43 * t52 + t49 * t62)) * MDP(18);];
taug = t1;
