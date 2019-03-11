% Calculate Gravitation load on the joints for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:41:18
% EndTime: 2019-03-08 19:41:20
% DurationCPUTime: 0.66s
% Computational Cost: add. (358->90), mult. (581->147), div. (0->0), fcn. (672->14), ass. (0->44)
t113 = MDP(15) - MDP(18);
t86 = sin(qJ(2));
t87 = cos(qJ(2));
t100 = cos(pkin(10));
t101 = cos(pkin(6));
t94 = t101 * t100;
t99 = sin(pkin(10));
t64 = t86 * t99 - t87 * t94;
t93 = t101 * t99;
t66 = t100 * t86 + t87 * t93;
t112 = -g(1) * t66 - g(2) * t64;
t82 = sin(pkin(6));
t109 = g(3) * t82;
t78 = pkin(12) + qJ(6);
t74 = sin(t78);
t79 = pkin(11) + qJ(4);
t77 = cos(t79);
t108 = t74 * t77;
t76 = cos(t78);
t107 = t76 * t77;
t80 = sin(pkin(12));
t106 = t77 * t80;
t83 = cos(pkin(12));
t105 = t77 * t83;
t104 = t77 * t87;
t103 = t82 * t86;
t102 = t82 * t87;
t98 = MDP(19) + MDP(8);
t97 = t82 * t100;
t96 = t82 * t99;
t65 = t86 * t94 + t87 * t99;
t67 = t100 * t87 - t86 * t93;
t95 = g(1) * t67 + g(2) * t65;
t75 = sin(t79);
t58 = t65 * t75 + t77 * t97;
t60 = t67 * t75 - t77 * t96;
t62 = -t101 * t77 + t103 * t75;
t91 = g(1) * t60 + g(2) * t58 + g(3) * t62;
t89 = g(3) * t102 + t112;
t84 = cos(pkin(11));
t63 = t101 * t75 + t103 * t77;
t61 = t67 * t77 + t75 * t96;
t59 = t65 * t77 - t75 * t97;
t1 = [(-MDP(1) - t98) * g(3); (-g(1) * (-t66 * pkin(2) + t67 * qJ(3)) - g(2) * (-t64 * pkin(2) + t65 * qJ(3)) - (pkin(2) * t87 + qJ(3) * t86) * t109) * MDP(8) + (-g(1) * (-t105 * t66 + t67 * t80) - g(2) * (-t105 * t64 + t65 * t80) - (t104 * t83 + t80 * t86) * t109) * MDP(16) + (-g(1) * (t106 * t66 + t67 * t83) - g(2) * (t106 * t64 + t65 * t83) - (-t104 * t80 + t83 * t86) * t109) * MDP(17) + ((t109 * t86 + t95) * (-pkin(8) - qJ(3)) + (-t87 * t109 - t112) * (pkin(3) * t84 + pkin(4) * t77 + qJ(5) * t75 + pkin(2))) * MDP(19) + (-g(1) * (-t107 * t66 + t67 * t74) - g(2) * (-t107 * t64 + t65 * t74) - (t104 * t76 + t74 * t86) * t109) * MDP(25) + (-g(1) * (t108 * t66 + t67 * t76) - g(2) * (t108 * t64 + t65 * t76) - (-t104 * t74 + t76 * t86) * t109) * MDP(26) + (MDP(4) - MDP(7)) * (g(3) * t103 + t95) + (-t77 * MDP(14) - t84 * MDP(5) + MDP(6) * sin(pkin(11)) + t113 * t75 - MDP(3)) * t89; t98 * t89; (-g(1) * (-pkin(4) * t60 + qJ(5) * t61) - g(2) * (-pkin(4) * t58 + qJ(5) * t59) - g(3) * (-pkin(4) * t62 + qJ(5) * t63)) * MDP(19) + t113 * (g(1) * t61 + g(2) * t59 + g(3) * t63) + (MDP(16) * t83 - MDP(17) * t80 + MDP(25) * t76 - MDP(26) * t74 + MDP(14)) * t91; -t91 * MDP(19); (-g(1) * (-t61 * t74 + t66 * t76) - g(2) * (-t59 * t74 + t64 * t76) - g(3) * (-t102 * t76 - t63 * t74)) * MDP(25) + (-g(1) * (-t61 * t76 - t66 * t74) - g(2) * (-t59 * t76 - t64 * t74) - g(3) * (t102 * t74 - t63 * t76)) * MDP(26);];
taug  = t1;
