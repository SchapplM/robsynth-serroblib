% Calculate Gravitation load on the joints for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:11:07
% EndTime: 2019-03-09 00:11:09
% DurationCPUTime: 0.53s
% Computational Cost: add. (366->103), mult. (730->173), div. (0->0), fcn. (861->12), ass. (0->48)
t115 = MDP(11) - MDP(26);
t100 = sin(pkin(11));
t87 = sin(qJ(2));
t90 = cos(qJ(2));
t101 = cos(pkin(11));
t102 = cos(pkin(6));
t97 = t102 * t101;
t64 = t100 * t87 - t90 * t97;
t114 = g(2) * t64;
t65 = t100 * t90 + t87 * t97;
t113 = g(2) * t65;
t84 = sin(pkin(6));
t112 = g(3) * t84;
t83 = qJ(4) + qJ(5);
t79 = sin(t83);
t85 = sin(qJ(4));
t71 = pkin(4) * t85 + pkin(5) * t79;
t111 = pkin(8) + t71;
t89 = cos(qJ(3));
t110 = t79 * t89;
t80 = cos(t83);
t109 = t80 * t89;
t108 = t84 * t87;
t107 = t84 * t90;
t106 = t85 * t89;
t88 = cos(qJ(4));
t105 = t88 * t89;
t104 = t89 * t90;
t86 = sin(qJ(3));
t99 = t84 * t101;
t61 = t65 * t89 - t86 * t99;
t96 = t102 * t100;
t67 = t101 * t90 - t87 * t96;
t98 = t84 * t100;
t63 = t67 * t89 + t86 * t98;
t66 = t101 * t87 + t90 * t96;
t69 = t102 * t86 + t89 * t108;
t91 = -g(1) * (-t63 * t79 + t66 * t80) - g(2) * (-t61 * t79 + t64 * t80) - g(3) * (-t80 * t107 - t69 * t79);
t103 = t91 * MDP(24) + (-g(1) * (-t63 * t80 - t66 * t79) - g(2) * (-t61 * t80 - t64 * t79) - g(3) * (t79 * t107 - t69 * t80)) * MDP(25);
t72 = t88 * pkin(4) + pkin(5) * t80;
t70 = pkin(3) + t72;
t82 = -qJ(6) - pkin(10) - pkin(9);
t95 = t70 * t89 - t82 * t86 + pkin(2);
t60 = t65 * t86 + t89 * t99;
t62 = t67 * t86 - t89 * t98;
t68 = -t102 * t89 + t86 * t108;
t94 = g(1) * t62 + g(2) * t60 + g(3) * t68;
t1 = [(-MDP(1) - MDP(27)) * g(3); (g(1) * t67 + g(3) * t108 + t113) * MDP(4) + (-g(1) * (-t66 * t105 + t67 * t85) - g(2) * (-t64 * t105 + t65 * t85) - (t88 * t104 + t85 * t87) * t112) * MDP(17) + (-g(1) * (t66 * t106 + t67 * t88) - g(2) * (t64 * t106 + t65 * t88) - (-t85 * t104 + t87 * t88) * t112) * MDP(18) + (-g(1) * (-t66 * t109 + t67 * t79) - g(2) * (-t64 * t109 + t65 * t79) - (t80 * t104 + t79 * t87) * t112) * MDP(24) + (-g(1) * (t66 * t110 + t67 * t80) - g(2) * (t64 * t110 + t65 * t80) - (-t79 * t104 + t80 * t87) * t112) * MDP(25) + (-g(1) * (t111 * t67 - t95 * t66) - t111 * t113 + t95 * t114 - (t111 * t87 + t95 * t90) * t112) * MDP(27) + (-t89 * MDP(10) + t115 * t86 - MDP(3)) * (-g(1) * t66 + g(3) * t107 - t114); (-g(1) * (-t62 * t70 - t63 * t82) - g(2) * (-t60 * t70 - t61 * t82) - g(3) * (-t68 * t70 - t69 * t82)) * MDP(27) + t115 * (g(1) * t63 + g(2) * t61 + g(3) * t69) + (MDP(17) * t88 - MDP(18) * t85 + MDP(24) * t80 - MDP(25) * t79 + MDP(10)) * t94; (-g(1) * (-t63 * t85 + t66 * t88) - g(2) * (-t61 * t85 + t64 * t88) - g(3) * (-t88 * t107 - t69 * t85)) * MDP(17) + (-g(1) * (-t63 * t88 - t66 * t85) - g(2) * (-t61 * t88 - t64 * t85) - g(3) * (t85 * t107 - t69 * t88)) * MDP(18) + (-g(1) * (-t63 * t71 + t66 * t72) - g(2) * (-t61 * t71 + t64 * t72) - g(3) * (-t72 * t107 - t69 * t71)) * MDP(27) + t103; MDP(27) * pkin(5) * t91 + t103; -t94 * MDP(27);];
taug  = t1;
