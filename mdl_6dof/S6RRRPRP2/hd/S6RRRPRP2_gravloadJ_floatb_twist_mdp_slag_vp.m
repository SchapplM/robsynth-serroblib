% Calculate Gravitation load on the joints for
% S6RRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:37:19
% EndTime: 2019-03-09 16:37:20
% DurationCPUTime: 0.40s
% Computational Cost: add. (409->86), mult. (415->114), div. (0->0), fcn. (380->10), ass. (0->47)
t102 = sin(qJ(5));
t105 = cos(qJ(5));
t133 = pkin(5) * t105 + qJ(6) * t102;
t130 = MDP(25) + MDP(27);
t129 = MDP(26) - MDP(29);
t104 = sin(qJ(1));
t107 = cos(qJ(1));
t87 = g(1) * t107 + g(2) * t104;
t101 = qJ(2) + qJ(3);
t95 = pkin(10) + t101;
t90 = sin(t95);
t132 = t87 * t90;
t88 = t90 * pkin(9);
t91 = cos(t95);
t131 = -t91 * pkin(4) - t88;
t96 = sin(t101);
t128 = pkin(3) * t96;
t127 = g(3) * t90;
t97 = cos(t101);
t92 = pkin(3) * t97;
t106 = cos(qJ(2));
t98 = t106 * pkin(2);
t123 = t92 + t98;
t122 = t107 * t91;
t100 = -qJ(4) - pkin(8) - pkin(7);
t120 = t100 * t107;
t119 = t104 * t102;
t118 = t104 * t105;
t117 = t105 * t107;
t116 = t107 * t102;
t115 = t133 * t91 - t131 + t92;
t86 = g(1) * t104 - g(2) * t107;
t109 = -g(3) * t97 + t87 * t96;
t112 = (-t87 * t91 - t127) * MDP(28) + t109 * MDP(16) + (g(3) * t96 + t87 * t97) * MDP(17) + (-t129 * t102 + t130 * t105) * (-g(3) * t91 + t132);
t75 = t91 * t119 + t117;
t77 = t91 * t116 - t118;
t64 = g(1) * t77 + g(2) * t75 + t102 * t127;
t108 = (pkin(4) + t133) * t132;
t103 = sin(qJ(2));
t85 = pkin(9) * t122;
t83 = t104 * t91 * pkin(9);
t81 = -pkin(2) * t103 - t128;
t80 = pkin(1) + t123;
t79 = t107 * t80;
t78 = t91 * t117 + t119;
t76 = t91 * t118 - t116;
t1 = [(-g(1) * (-t104 * t80 - t120) - g(2) * (-t104 * t100 + t79)) * MDP(19) + (-g(1) * (-t76 * pkin(5) - t75 * qJ(6) - t120) - g(2) * (pkin(4) * t122 + t78 * pkin(5) + t77 * qJ(6) + t107 * t88 + t79) + (-g(1) * (-t80 + t131) + g(2) * t100) * t104) * MDP(30) + (MDP(3) - MDP(18)) * t87 + t130 * (g(1) * t76 - g(2) * t78) - t129 * (g(1) * t75 - g(2) * t77) + (-t103 * MDP(10) + MDP(16) * t97 - MDP(17) * t96 + t90 * MDP(28) + t106 * MDP(9) + MDP(2)) * t86; (-g(3) * t106 + t87 * t103) * MDP(9) + (g(3) * t103 + t87 * t106) * MDP(10) + (-g(3) * t123 - t87 * t81) * MDP(19) + (-g(1) * (t107 * t81 + t85) - g(2) * (t104 * t81 + t83) - g(3) * (t98 + t115) + t108) * MDP(30) + t112; t109 * pkin(3) * MDP(19) + (-g(1) * (-t107 * t128 + t85) - g(2) * (-t104 * t128 + t83) - g(3) * t115 + t108) * MDP(30) + t112; (-MDP(19) - MDP(30)) * t86; (-g(1) * (-pkin(5) * t77 + qJ(6) * t78) - g(2) * (-pkin(5) * t75 + qJ(6) * t76) - (-pkin(5) * t102 + qJ(6) * t105) * t127) * MDP(30) + t130 * t64 + t129 * (g(1) * t78 + g(2) * t76 + t105 * t127); -t64 * MDP(30);];
taug  = t1;
