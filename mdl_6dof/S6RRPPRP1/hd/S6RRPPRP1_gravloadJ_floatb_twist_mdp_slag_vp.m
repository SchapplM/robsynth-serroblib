% Calculate Gravitation load on the joints for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPPRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:27:48
% EndTime: 2019-03-09 08:27:50
% DurationCPUTime: 0.70s
% Computational Cost: add. (362->94), mult. (392->124), div. (0->0), fcn. (364->10), ass. (0->46)
t121 = MDP(15) + MDP(25);
t120 = MDP(22) + MDP(24);
t119 = MDP(23) - MDP(26);
t90 = sin(qJ(1));
t92 = cos(qJ(1));
t71 = g(1) * t92 + g(2) * t90;
t84 = qJ(2) + pkin(9);
t78 = sin(t84);
t80 = cos(t84);
t96 = -g(3) * t80 + t71 * t78;
t89 = sin(qJ(2));
t118 = pkin(2) * t89;
t115 = g(3) * t78;
t88 = -pkin(8) - qJ(4);
t113 = t78 * t88;
t83 = pkin(10) + qJ(5);
t79 = cos(t83);
t112 = t79 * t92;
t85 = sin(pkin(10));
t111 = t85 * t92;
t86 = cos(pkin(10));
t110 = t86 * t92;
t77 = sin(t83);
t109 = t90 * t77;
t108 = t90 * t79;
t107 = t90 * t85;
t106 = t90 * t86;
t105 = t92 * t77;
t104 = MDP(16) + MDP(27);
t87 = -qJ(3) - pkin(7);
t103 = pkin(4) * t85 - t87;
t70 = g(1) * t90 - g(2) * t92;
t101 = pkin(3) * t80 + qJ(4) * t78;
t75 = pkin(4) * t86 + pkin(3);
t100 = t75 * t80 - t113;
t99 = pkin(5) * t79 + qJ(6) * t77 + t75;
t65 = t80 * t109 + t112;
t67 = t80 * t105 - t108;
t59 = g(1) * t67 + g(2) * t65 + t77 * t115;
t91 = cos(qJ(2));
t81 = t91 * pkin(2);
t76 = t81 + pkin(1);
t72 = t92 * t76;
t68 = t80 * t112 + t109;
t66 = t80 * t108 - t105;
t1 = [(-g(1) * (-t90 * t76 - t92 * t87) - g(2) * (-t90 * t87 + t72)) * MDP(12) + (-g(1) * (-t80 * t106 + t111) - g(2) * (t80 * t110 + t107)) * MDP(13) + (-g(1) * (t80 * t107 + t110) - g(2) * (-t80 * t111 + t106)) * MDP(14) + (-g(2) * t72 + (g(1) * t87 - g(2) * t101) * t92 + (-g(1) * (-t101 - t76) + g(2) * t87) * t90) * MDP(16) + (-g(1) * (-t66 * pkin(5) - t65 * qJ(6)) - g(2) * (t68 * pkin(5) + t67 * qJ(6) + t72) + (-g(1) * t103 - g(2) * t100) * t92 + (-g(1) * (-t100 - t76) - g(2) * t103) * t90) * MDP(27) + (MDP(3) - MDP(11)) * t71 + t120 * (g(1) * t66 - g(2) * t68) - t119 * (g(1) * t65 - g(2) * t67) + (-t89 * MDP(10) + t91 * MDP(9) + t121 * t78 + MDP(2)) * t70; (g(3) * t89 + t71 * t91) * MDP(10) + (-g(3) * (t101 + t81) + t71 * (pkin(3) * t78 - qJ(4) * t80 + t118)) * MDP(16) + (-g(3) * (t99 * t80 - t113 + t81) + t71 * (t99 * t78 + t80 * t88 + t118)) * MDP(27) + (pkin(2) * MDP(12) + MDP(9)) * (-g(3) * t91 + t71 * t89) + t121 * (-t71 * t80 - t115) + (t86 * MDP(13) - t85 * MDP(14) - t119 * t77 + t120 * t79) * t96; (-MDP(12) - t104) * t70; -t104 * t96; (-g(1) * (-pkin(5) * t67 + qJ(6) * t68) - g(2) * (-pkin(5) * t65 + qJ(6) * t66) - (-pkin(5) * t77 + qJ(6) * t79) * t115) * MDP(27) + t119 * (g(1) * t68 + g(2) * t66 + t79 * t115) + t120 * t59; -t59 * MDP(27);];
taug  = t1;
