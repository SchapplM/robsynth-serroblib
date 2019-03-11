% Calculate Gravitation load on the joints for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:25:06
% EndTime: 2019-03-08 20:25:08
% DurationCPUTime: 0.42s
% Computational Cost: add. (342->74), mult. (708->138), div. (0->0), fcn. (887->14), ass. (0->39)
t105 = cos(pkin(12));
t85 = sin(pkin(12));
t92 = sin(qJ(2));
t95 = cos(qJ(2));
t100 = t95 * t105 - t92 * t85;
t87 = sin(pkin(6));
t115 = g(3) * t87;
t84 = qJ(4) + qJ(5);
t83 = cos(t84);
t90 = sin(qJ(6));
t114 = t83 * t90;
t93 = cos(qJ(6));
t113 = t83 * t93;
t86 = sin(pkin(11));
t112 = t86 * t87;
t88 = cos(pkin(11));
t111 = t87 * t88;
t91 = sin(qJ(4));
t110 = t87 * t91;
t94 = cos(qJ(4));
t109 = t87 * t94;
t89 = cos(pkin(6));
t108 = t89 * t92;
t107 = t89 * t95;
t79 = -t92 * t105 - t95 * t85;
t77 = t79 * t89;
t101 = t100 * t88 + t77 * t86;
t102 = t100 * t86 - t77 * t88;
t82 = sin(t84);
t63 = t102 * t83 - t82 * t111;
t65 = t101 * t83 + t82 * t112;
t76 = t79 * t87;
t73 = -t76 * t83 + t82 * t89;
t103 = (g(1) * t65 + g(2) * t63 + g(3) * t73) * MDP(19) + (-t93 * MDP(25) + t90 * MDP(26) - MDP(18)) * (g(1) * (-t101 * t82 + t83 * t112) + g(2) * (-t102 * t82 - t83 * t111) + g(3) * (t76 * t82 + t83 * t89));
t97 = t89 * t100;
t75 = t100 * t87;
t70 = t79 * t88 - t86 * t97;
t67 = t86 * t79 + t88 * t97;
t1 = [(-MDP(1) - MDP(5)) * g(3); (-g(1) * (t86 * t108 - t88 * t95) - g(2) * (-t88 * t108 - t86 * t95) + t92 * t115) * MDP(4) + (-g(1) * (t101 * t90 + t70 * t113) - g(2) * (t102 * t90 + t67 * t113) - g(3) * (t75 * t113 - t76 * t90)) * MDP(25) + (-g(1) * (t101 * t93 - t70 * t114) - g(2) * (t102 * t93 - t67 * t114) - g(3) * (-t75 * t114 - t76 * t93)) * MDP(26) + (pkin(2) * MDP(5) + MDP(3)) * (-g(1) * (-t86 * t107 - t88 * t92) - g(2) * (t88 * t107 - t86 * t92) - t95 * t115) + (-MDP(11) * t94 + MDP(12) * t91 - t83 * MDP(18) + MDP(19) * t82) * (g(1) * t70 + g(2) * t67 + g(3) * t75); (-g(3) * t89 + (-g(1) * t86 + g(2) * t88) * t87) * MDP(5); (-g(1) * (-t101 * t91 + t86 * t109) - g(2) * (-t102 * t91 - t88 * t109) - g(3) * (t76 * t91 + t89 * t94)) * MDP(11) + (-g(1) * (-t101 * t94 - t86 * t110) - g(2) * (-t102 * t94 + t88 * t110) - g(3) * (t76 * t94 - t89 * t91)) * MDP(12) + t103; t103; (-g(1) * (-t65 * t90 - t70 * t93) - g(2) * (-t63 * t90 - t67 * t93) - g(3) * (-t73 * t90 - t75 * t93)) * MDP(25) + (-g(1) * (-t65 * t93 + t70 * t90) - g(2) * (-t63 * t93 + t67 * t90) - g(3) * (-t73 * t93 + t75 * t90)) * MDP(26);];
taug  = t1;
