% Calculate Gravitation load on the joints for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPPRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RPPPRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:32:12
% EndTime: 2019-03-09 01:32:12
% DurationCPUTime: 0.17s
% Computational Cost: add. (147->46), mult. (142->64), div. (0->0), fcn. (118->10), ass. (0->27)
t50 = pkin(10) + qJ(5);
t47 = cos(t50);
t70 = MDP(18) * t47;
t54 = sin(qJ(6));
t56 = cos(qJ(6));
t69 = t56 * MDP(24) - t54 * MDP(25) + MDP(17);
t68 = g(3) * t47;
t51 = qJ(1) + pkin(9);
t46 = sin(t51);
t67 = t46 * t54;
t66 = t46 * t56;
t48 = cos(t51);
t65 = t48 * t54;
t64 = t48 * t56;
t61 = -MDP(11) - MDP(7);
t57 = cos(qJ(1));
t60 = t57 * pkin(1) + t48 * pkin(2) + t46 * qJ(3);
t55 = sin(qJ(1));
t59 = -t55 * pkin(1) + t48 * qJ(3);
t40 = -g(1) * t48 - g(2) * t46;
t39 = g(1) * t46 - g(2) * t48;
t45 = sin(t50);
t38 = t45 * t64 - t67;
t37 = t45 * t65 + t66;
t36 = t45 * t66 + t65;
t35 = -t45 * t67 + t64;
t1 = [(g(1) * t57 + g(2) * t55) * MDP(3) + (-g(1) * (-t46 * pkin(2) + t59) - g(2) * t60) * MDP(7) + (-g(1) * ((-pkin(2) - qJ(4)) * t46 + t59) - g(2) * (t48 * qJ(4) + t60)) * MDP(11) + (-g(1) * t38 - g(2) * t36) * MDP(24) + (g(1) * t37 - g(2) * t35) * MDP(25) + (-MDP(5) + MDP(10)) * t39 + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t55 - g(2) * t57) + (MDP(17) * t45 + t70 + MDP(8) * sin(pkin(10)) + MDP(9) * cos(pkin(10)) + MDP(6)) * t40; (-MDP(4) + t61) * g(3); t61 * t39; t40 * MDP(11); (t69 * t45 + t70) * g(3) + (MDP(18) * t45 - t69 * t47) * t39; (-g(1) * t35 - g(2) * t37 + t54 * t68) * MDP(24) + (g(1) * t36 - g(2) * t38 + t56 * t68) * MDP(25);];
taug  = t1;
