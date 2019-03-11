% Calculate Gravitation load on the joints for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPPRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPPRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:16:15
% EndTime: 2019-03-08 19:16:16
% DurationCPUTime: 0.52s
% Computational Cost: add. (289->78), mult. (627->139), div. (0->0), fcn. (772->14), ass. (0->41)
t102 = cos(pkin(11));
t83 = sin(pkin(11));
t90 = sin(qJ(2));
t92 = cos(qJ(2));
t97 = t102 * t92 - t83 * t90;
t81 = pkin(12) + qJ(5);
t80 = cos(t81);
t89 = sin(qJ(6));
t111 = t80 * t89;
t91 = cos(qJ(6));
t110 = t80 * t91;
t84 = sin(pkin(10));
t85 = sin(pkin(6));
t109 = t84 * t85;
t108 = t84 * t90;
t87 = cos(pkin(10));
t107 = t85 * t87;
t106 = t85 * t92;
t88 = cos(pkin(6));
t105 = t88 * t90;
t104 = t88 * t92;
t101 = MDP(5) + MDP(9);
t100 = t87 * t104;
t75 = -t102 * t90 - t83 * t92;
t73 = t75 * t88;
t62 = -t73 * t87 + t84 * t97;
t63 = -t73 * t84 - t87 * t97;
t98 = -t104 * t84 - t87 * t90;
t94 = t88 * t97;
t61 = t84 * t75 + t87 * t94;
t64 = t75 * t87 - t84 * t94;
t71 = t97 * t85;
t95 = g(1) * t64 + g(2) * t61 + g(3) * t71;
t93 = -g(1) * t98 - g(3) * t106;
t79 = sin(t81);
t76 = pkin(2) * t100;
t72 = t75 * t85;
t67 = -t72 * t80 + t79 * t88;
t59 = t109 * t79 - t63 * t80;
t57 = -t107 * t79 + t62 * t80;
t1 = [(-MDP(1) - t101) * g(3); (-g(2) * (t100 - t108) + t93) * MDP(3) + (-g(1) * (t105 * t84 - t87 * t92) - g(2) * (-t105 * t87 - t84 * t92) + g(3) * t85 * t90) * MDP(4) + (-g(2) * t76 + (g(2) * t108 + t93) * pkin(2)) * MDP(5) + (g(1) * t63 - g(2) * t62 + g(3) * t72) * MDP(8) + (-g(1) * (pkin(2) * t98 + t64 * pkin(3) - t63 * qJ(4)) - g(2) * (-pkin(2) * t108 + pkin(3) * t61 + qJ(4) * t62 + t76) - g(3) * (pkin(2) * t106 + pkin(3) * t71 - qJ(4) * t72)) * MDP(9) + (-g(1) * (t110 * t64 - t63 * t89) - g(2) * (t110 * t61 + t62 * t89) - g(3) * (t110 * t71 - t72 * t89)) * MDP(22) + (-g(1) * (-t111 * t64 - t63 * t91) - g(2) * (-t111 * t61 + t62 * t91) - g(3) * (-t111 * t71 - t72 * t91)) * MDP(23) + (-t80 * MDP(15) + MDP(16) * t79 - MDP(6) * cos(pkin(12)) + MDP(7) * sin(pkin(12))) * t95; t101 * (-g(3) * t88 + (-g(1) * t84 + g(2) * t87) * t85); t95 * MDP(9); (g(1) * t59 + g(2) * t57 + g(3) * t67) * MDP(16) + (-MDP(22) * t91 + MDP(23) * t89 - MDP(15)) * (g(1) * (t109 * t80 + t63 * t79) + g(2) * (-t107 * t80 - t62 * t79) + g(3) * (t72 * t79 + t80 * t88)); (-g(1) * (-t59 * t89 - t64 * t91) - g(2) * (-t57 * t89 - t61 * t91) - g(3) * (-t67 * t89 - t71 * t91)) * MDP(22) + (-g(1) * (-t59 * t91 + t64 * t89) - g(2) * (-t57 * t91 + t61 * t89) - g(3) * (-t67 * t91 + t71 * t89)) * MDP(23);];
taug  = t1;
