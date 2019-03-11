% Calculate Gravitation load on the joints for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPPRRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:23:49
% EndTime: 2019-03-09 02:23:49
% DurationCPUTime: 0.21s
% Computational Cost: add. (182->54), mult. (191->84), div. (0->0), fcn. (180->10), ass. (0->30)
t56 = cos(qJ(4));
t72 = MDP(14) * t56;
t51 = qJ(5) + qJ(6);
t48 = sin(t51);
t49 = cos(t51);
t52 = sin(qJ(5));
t55 = cos(qJ(5));
t71 = t55 * MDP(20) - t52 * MDP(21) + t49 * MDP(27) - t48 * MDP(28) + MDP(13);
t70 = g(3) * t56;
t53 = sin(qJ(4));
t69 = t48 * t53;
t68 = t49 * t53;
t67 = t52 * t53;
t66 = t53 * t55;
t50 = qJ(1) + pkin(10);
t46 = sin(t50);
t47 = cos(t50);
t37 = -t46 * t69 + t47 * t49;
t38 = t46 * t68 + t47 * t48;
t39 = t46 * t49 + t47 * t69;
t40 = -t46 * t48 + t47 * t68;
t65 = (-g(1) * t37 - g(2) * t39 + t48 * t70) * MDP(27) + (g(1) * t38 - g(2) * t40 + t49 * t70) * MDP(28);
t59 = g(1) * t46 - g(2) * t47;
t57 = cos(qJ(1));
t54 = sin(qJ(1));
t44 = -t46 * t52 + t47 * t66;
t43 = t46 * t55 + t47 * t67;
t42 = t46 * t66 + t47 * t52;
t41 = -t46 * t67 + t47 * t55;
t1 = [(g(1) * t57 + g(2) * t54) * MDP(3) - t59 * MDP(5) + (-g(1) * (-t54 * pkin(1) - t46 * pkin(2) + t47 * qJ(3)) - g(2) * (t57 * pkin(1) + t47 * pkin(2) + t46 * qJ(3))) * MDP(7) + (-g(1) * t44 - g(2) * t42) * MDP(20) + (g(1) * t43 - g(2) * t41) * MDP(21) + (-g(1) * t40 - g(2) * t38) * MDP(27) + (g(1) * t39 - g(2) * t37) * MDP(28) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t54 - g(2) * t57) + (MDP(13) * t53 + MDP(6) + t72) * (-g(1) * t47 - g(2) * t46); (-MDP(4) - MDP(7)) * g(3); -t59 * MDP(7); (t71 * t53 + t72) * g(3) + (MDP(14) * t53 - t71 * t56) * t59; (-g(1) * t41 - g(2) * t43 + t52 * t70) * MDP(20) + (g(1) * t42 - g(2) * t44 + t55 * t70) * MDP(21) + t65; t65;];
taug  = t1;
