% Calculate Gravitation load on the joints for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:59:14
% EndTime: 2019-03-09 21:59:15
% DurationCPUTime: 0.35s
% Computational Cost: add. (441->76), mult. (377->106), div. (0->0), fcn. (331->12), ass. (0->45)
t93 = qJ(2) + qJ(3);
t89 = qJ(4) + t93;
t83 = sin(t89);
t84 = cos(t89);
t110 = t84 * pkin(4) + t83 * qJ(5);
t88 = cos(t93);
t108 = pkin(3) * t88 + t110;
t98 = cos(qJ(2));
t124 = t98 * pkin(2) + t108;
t123 = MDP(24) - MDP(27);
t121 = pkin(4) * t83;
t87 = sin(t93);
t106 = -pkin(3) * t87 - t121;
t97 = sin(qJ(1));
t99 = cos(qJ(1));
t105 = g(1) * t99 + g(2) * t97;
t67 = -g(3) * t84 + t105 * t83;
t120 = g(3) * t83;
t91 = pkin(11) + qJ(6);
t85 = sin(t91);
t118 = t97 * t85;
t86 = cos(t91);
t117 = t97 * t86;
t94 = sin(pkin(11));
t116 = t97 * t94;
t95 = cos(pkin(11));
t115 = t97 * t95;
t114 = t99 * t85;
t113 = t99 * t86;
t112 = t99 * t94;
t111 = t99 * t95;
t109 = qJ(5) * t84;
t96 = sin(qJ(2));
t107 = -t96 * pkin(2) + t106;
t103 = pkin(1) + t124;
t102 = t123 * (t105 * t84 + t120) + (MDP(25) * t95 - MDP(26) * t94 + MDP(34) * t86 - MDP(35) * t85 + MDP(23)) * t67;
t101 = (-g(3) * t88 + t105 * t87) * MDP(16) + (g(3) * t87 + t105 * t88) * MDP(17) + t102;
t92 = -pkin(9) - pkin(8) - pkin(7);
t79 = t99 * t109;
t78 = t97 * t109;
t74 = t84 * t113 + t118;
t73 = -t84 * t114 + t117;
t72 = -t84 * t117 + t114;
t71 = t84 * t118 + t113;
t1 = [t105 * MDP(3) + (-g(1) * (-t84 * t115 + t112) - g(2) * (t84 * t111 + t116)) * MDP(25) + (-g(1) * (t84 * t116 + t111) - g(2) * (-t84 * t112 + t115)) * MDP(26) + ((g(1) * t92 - g(2) * t103) * t99 + (g(1) * t103 + g(2) * t92) * t97) * MDP(28) + (-g(1) * t72 - g(2) * t74) * MDP(34) + (-g(1) * t71 - g(2) * t73) * MDP(35) + (-t96 * MDP(10) + MDP(16) * t88 - MDP(17) * t87 + t84 * MDP(23) + t98 * MDP(9) - t123 * t83 + MDP(2)) * (g(1) * t97 - g(2) * t99); (-g(3) * t98 + t105 * t96) * MDP(9) + (g(3) * t96 + t105 * t98) * MDP(10) + (-g(1) * (t107 * t99 + t79) - g(2) * (t107 * t97 + t78) - g(3) * t124) * MDP(28) + t101; (-g(1) * (t106 * t99 + t79) - g(2) * (t106 * t97 + t78) - g(3) * t108) * MDP(28) + t101; (-g(1) * (-t99 * t121 + t79) - g(2) * (-t97 * t121 + t78) - g(3) * t110) * MDP(28) + t102; -t67 * MDP(28); (-g(1) * t73 + g(2) * t71 + t85 * t120) * MDP(34) + (g(1) * t74 - g(2) * t72 + t86 * t120) * MDP(35);];
taug  = t1;
