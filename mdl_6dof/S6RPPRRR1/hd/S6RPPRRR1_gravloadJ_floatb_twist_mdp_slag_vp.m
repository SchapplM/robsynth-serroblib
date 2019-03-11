% Calculate Gravitation load on the joints for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:18:39
% EndTime: 2019-03-09 02:18:39
% DurationCPUTime: 0.16s
% Computational Cost: add. (214->46), mult. (170->66), div. (0->0), fcn. (146->12), ass. (0->26)
t48 = pkin(11) + qJ(4);
t47 = qJ(5) + t48;
t41 = sin(t47);
t66 = g(3) * t41;
t49 = qJ(1) + pkin(10);
t44 = sin(t49);
t52 = sin(qJ(6));
t64 = t44 * t52;
t54 = cos(qJ(6));
t63 = t44 * t54;
t46 = cos(t49);
t62 = t46 * t52;
t61 = t46 * t54;
t42 = cos(t47);
t59 = g(1) * t46 + g(2) * t44;
t60 = (t59 * t42 + t66) * MDP(22) + (t54 * MDP(28) - t52 * MDP(29) + MDP(21)) * (-g(3) * t42 + t59 * t41);
t58 = g(1) * t44 - g(2) * t46;
t55 = cos(qJ(1));
t53 = sin(qJ(1));
t45 = cos(t48);
t43 = sin(t48);
t40 = t42 * t61 + t64;
t39 = -t42 * t62 + t63;
t38 = -t42 * t63 + t62;
t37 = t42 * t64 + t61;
t1 = [(g(1) * t55 + g(2) * t53) * MDP(3) - t59 * MDP(7) + (-g(1) * (-t53 * pkin(1) - t44 * pkin(2) + t46 * qJ(3)) - g(2) * (t55 * pkin(1) + t46 * pkin(2) + t44 * qJ(3))) * MDP(8) + (-g(1) * t38 - g(2) * t40) * MDP(28) + (-g(1) * t37 - g(2) * t39) * MDP(29) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t53 - g(2) * t55) + (t45 * MDP(14) - t43 * MDP(15) + MDP(21) * t42 - MDP(22) * t41 + MDP(5) * cos(pkin(11)) - MDP(6) * sin(pkin(11))) * t58; (-MDP(4) - MDP(8)) * g(3); -t58 * MDP(8); (-g(3) * t45 + t59 * t43) * MDP(14) + (g(3) * t43 + t59 * t45) * MDP(15) + t60; t60; (-g(1) * t39 + g(2) * t37 + t52 * t66) * MDP(28) + (g(1) * t40 - g(2) * t38 + t54 * t66) * MDP(29);];
taug  = t1;
