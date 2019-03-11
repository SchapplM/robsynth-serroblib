% Calculate Gravitation load on the joints for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPPRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:48:04
% EndTime: 2019-03-09 08:48:05
% DurationCPUTime: 0.34s
% Computational Cost: add. (244->65), mult. (336->95), div. (0->0), fcn. (339->10), ass. (0->32)
t74 = sin(qJ(1));
t69 = qJ(2) + pkin(10);
t65 = sin(t69);
t66 = cos(t69);
t72 = sin(qJ(5));
t90 = cos(qJ(5));
t95 = t65 * t90 - t66 * t72;
t51 = t95 * t74;
t82 = t65 * t72 + t66 * t90;
t52 = t82 * t74;
t77 = cos(qJ(1));
t53 = t95 * t77;
t54 = t82 * t77;
t71 = sin(qJ(6));
t75 = cos(qJ(6));
t91 = g(3) * t95;
t97 = (MDP(29) * t75 - MDP(30) * t71 + MDP(22)) * (-g(1) * t53 - g(2) * t51 + g(3) * t82) + (g(1) * t54 + g(2) * t52 + t91) * MDP(23);
t60 = g(1) * t77 + g(2) * t74;
t59 = g(1) * t74 - g(2) * t77;
t85 = t66 * pkin(3) + t65 * qJ(4);
t84 = t52 * t75 + t77 * t71;
t83 = t52 * t71 - t77 * t75;
t76 = cos(qJ(2));
t73 = sin(qJ(2));
t70 = -qJ(3) - pkin(7);
t67 = t76 * pkin(2);
t64 = t67 + pkin(1);
t61 = t77 * t64;
t50 = -g(3) * t66 + t60 * t65;
t49 = t54 * t75 - t74 * t71;
t48 = -t54 * t71 - t74 * t75;
t1 = [(-g(1) * (-t74 * t64 - t77 * t70) - g(2) * (-t74 * t70 + t61)) * MDP(12) + (-g(2) * t61 + (g(1) * t70 - g(2) * t85) * t77 + (-g(1) * (-t64 - t85) + g(2) * t70) * t74) * MDP(16) + (g(1) * t52 - g(2) * t54) * MDP(22) + (g(1) * t51 - g(2) * t53) * MDP(23) + (g(1) * t84 - g(2) * t49) * MDP(29) + (-g(1) * t83 - g(2) * t48) * MDP(30) + (MDP(3) - MDP(11) - MDP(14)) * t60 + (-t73 * MDP(10) + MDP(13) * t66 + t65 * MDP(15) + t76 * MDP(9) + MDP(2)) * t59; (g(3) * t73 + t60 * t76) * MDP(10) + t50 * MDP(13) + (-g(3) * t65 - t60 * t66) * MDP(15) + (-g(3) * (t67 + t85) + t60 * (pkin(2) * t73 + pkin(3) * t65 - qJ(4) * t66)) * MDP(16) + (pkin(2) * MDP(12) + MDP(9)) * (-g(3) * t76 + t60 * t73) - t97; (-MDP(12) - MDP(16)) * t59; -t50 * MDP(16); t97; (-g(1) * t48 + g(2) * t83 + t71 * t91) * MDP(29) + (g(1) * t49 + g(2) * t84 + t75 * t91) * MDP(30);];
taug  = t1;
