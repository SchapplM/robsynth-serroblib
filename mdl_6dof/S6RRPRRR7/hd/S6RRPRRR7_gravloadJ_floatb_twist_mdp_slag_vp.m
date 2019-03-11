% Calculate Gravitation load on the joints for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:59:46
% EndTime: 2019-03-09 13:59:47
% DurationCPUTime: 0.34s
% Computational Cost: add. (229->65), mult. (476->99), div. (0->0), fcn. (519->10), ass. (0->35)
t74 = sin(qJ(4));
t75 = sin(qJ(2));
t78 = cos(qJ(2));
t97 = cos(qJ(4));
t105 = -t78 * t74 + t75 * t97;
t76 = sin(qJ(1));
t59 = t105 * t76;
t65 = t75 * t74 + t78 * t97;
t60 = t65 * t76;
t79 = cos(qJ(1));
t61 = t105 * t79;
t62 = t65 * t79;
t72 = qJ(5) + qJ(6);
t70 = sin(t72);
t71 = cos(t72);
t73 = sin(qJ(5));
t77 = cos(qJ(5));
t98 = g(3) * t105;
t106 = (MDP(27) * t77 - MDP(28) * t73 + MDP(34) * t71 - MDP(35) * t70 + MDP(20)) * (-g(1) * t61 - g(2) * t59 + g(3) * t65) + (g(1) * t62 + g(2) * t60 + t98) * MDP(21);
t67 = g(1) * t79 + g(2) * t76;
t103 = MDP(9) + MDP(11);
t102 = MDP(10) - MDP(13);
t53 = -t62 * t70 - t76 * t71;
t54 = t62 * t71 - t76 * t70;
t85 = t60 * t70 - t79 * t71;
t86 = t60 * t71 + t79 * t70;
t95 = (-g(1) * t53 + g(2) * t85 + t70 * t98) * MDP(34) + (g(1) * t54 + g(2) * t86 + t71 * t98) * MDP(35);
t88 = t78 * pkin(2) + t75 * qJ(3);
t84 = t60 * t77 + t79 * t73;
t83 = t60 * t73 - t79 * t77;
t82 = pkin(1) + t88;
t57 = -g(3) * t78 + t67 * t75;
t56 = t62 * t77 - t76 * t73;
t55 = -t62 * t73 - t76 * t77;
t1 = [((-g(1) * pkin(7) - g(2) * t82) * t79 + (-g(2) * pkin(7) + g(1) * t82) * t76) * MDP(14) + (g(1) * t60 - g(2) * t62) * MDP(20) + (g(1) * t59 - g(2) * t61) * MDP(21) + (g(1) * t84 - g(2) * t56) * MDP(27) + (-g(1) * t83 - g(2) * t55) * MDP(28) + (g(1) * t86 - g(2) * t54) * MDP(34) + (-g(1) * t85 - g(2) * t53) * MDP(35) + (MDP(3) - MDP(12)) * t67 + (-t102 * t75 + t103 * t78 + MDP(2)) * (g(1) * t76 - g(2) * t79); (-g(3) * t88 + t67 * (pkin(2) * t75 - qJ(3) * t78)) * MDP(14) + t102 * (g(3) * t75 + t67 * t78) + t103 * t57 - t106; -t57 * MDP(14); t106; (-g(1) * t55 + g(2) * t83 + t73 * t98) * MDP(27) + (g(1) * t56 + g(2) * t84 + t77 * t98) * MDP(28) + t95; t95;];
taug  = t1;
