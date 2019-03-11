% Calculate Gravitation load on the joints for
% S6RPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:55:17
% EndTime: 2019-03-09 06:55:17
% DurationCPUTime: 0.14s
% Computational Cost: add. (266->42), mult. (208->61), div. (0->0), fcn. (182->12), ass. (0->28)
t53 = qJ(3) + qJ(4);
t51 = qJ(5) + t53;
t45 = sin(t51);
t71 = g(3) * t45;
t52 = qJ(1) + pkin(11);
t47 = sin(t52);
t54 = sin(qJ(6));
t69 = t47 * t54;
t57 = cos(qJ(6));
t68 = t47 * t57;
t48 = cos(t52);
t67 = t48 * t54;
t66 = t48 * t57;
t46 = cos(t51);
t64 = g(1) * t48 + g(2) * t47;
t65 = (t46 * t64 + t71) * MDP(25) + (t57 * MDP(31) - t54 * MDP(32) + MDP(24)) * (-g(3) * t46 + t64 * t45);
t49 = sin(t53);
t50 = cos(t53);
t61 = (-g(3) * t50 + t49 * t64) * MDP(17) + (g(3) * t49 + t50 * t64) * MDP(18) + t65;
t59 = cos(qJ(1));
t58 = cos(qJ(3));
t56 = sin(qJ(1));
t55 = sin(qJ(3));
t44 = t46 * t66 + t69;
t43 = -t46 * t67 + t68;
t42 = -t46 * t68 + t67;
t41 = t46 * t69 + t66;
t1 = [(g(1) * t59 + g(2) * t56) * MDP(3) + (-g(1) * t42 - g(2) * t44) * MDP(31) + (-g(1) * t41 - g(2) * t43) * MDP(32) + (MDP(4) * pkin(1) + MDP(2)) * (g(1) * t56 - g(2) * t59) + (t58 * MDP(10) - t55 * MDP(11) + MDP(17) * t50 - MDP(18) * t49 + MDP(24) * t46 - MDP(25) * t45) * (g(1) * t47 - g(2) * t48); -g(3) * MDP(4); (-g(3) * t58 + t55 * t64) * MDP(10) + (g(3) * t55 + t58 * t64) * MDP(11) + t61; t61; t65; (-g(1) * t43 + g(2) * t41 + t54 * t71) * MDP(31) + (g(1) * t44 - g(2) * t42 + t57 * t71) * MDP(32);];
taug  = t1;
