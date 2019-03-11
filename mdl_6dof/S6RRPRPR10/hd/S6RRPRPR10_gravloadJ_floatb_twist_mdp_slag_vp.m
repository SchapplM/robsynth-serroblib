% Calculate Gravitation load on the joints for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRPR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:09:26
% EndTime: 2019-03-09 11:09:30
% DurationCPUTime: 1.08s
% Computational Cost: add. (448->117), mult. (792->183), div. (0->0), fcn. (916->12), ass. (0->50)
t147 = MDP(10) - MDP(13) - MDP(22);
t146 = MDP(20) - MDP(23);
t141 = MDP(21) - MDP(24);
t107 = sin(qJ(6));
t110 = cos(qJ(6));
t102 = pkin(11) + qJ(4);
t100 = cos(t102);
t104 = sin(pkin(6));
t138 = cos(qJ(1));
t125 = t104 * t138;
t108 = sin(qJ(2));
t109 = sin(qJ(1));
t111 = cos(qJ(2));
t131 = cos(pkin(6));
t117 = t131 * t138;
t88 = t108 * t117 + t109 * t111;
t99 = sin(t102);
t79 = t100 * t125 + t88 * t99;
t87 = t108 * t109 - t111 * t117;
t143 = t107 * t79 + t110 * t87;
t142 = -t107 * t87 + t110 * t79;
t137 = g(3) * t104;
t123 = t109 * t131;
t89 = t108 * t138 + t111 * t123;
t112 = -g(1) * t89 - g(2) * t87 + t111 * t137;
t135 = t107 * t99;
t133 = t110 * t99;
t129 = t104 * t109;
t132 = t138 * pkin(1) + pkin(8) * t129;
t130 = t104 * t108;
t128 = t107 * t111;
t127 = t110 * t111;
t103 = sin(pkin(11));
t126 = t103 * t129;
t124 = -t109 * pkin(1) + pkin(8) * t125;
t80 = t100 * t88 - t99 * t125;
t122 = t103 * t125;
t90 = -t108 * t123 + t111 * t138;
t118 = g(1) * t90 + g(2) * t88;
t83 = -t100 * t129 + t90 * t99;
t85 = -t100 * t131 + t130 * t99;
t115 = g(1) * t83 + g(2) * t79 + g(3) * t85;
t106 = -pkin(9) - qJ(3);
t105 = cos(pkin(11));
t98 = pkin(3) * t105 + pkin(2);
t86 = t100 * t130 + t131 * t99;
t84 = t100 * t90 + t129 * t99;
t75 = t107 * t83 + t110 * t89;
t74 = -t107 * t89 + t110 * t83;
t1 = [(g(1) * t109 - g(2) * t138) * MDP(2) + (g(1) * t138 + g(2) * t109) * MDP(3) + (g(1) * t88 - g(2) * t90) * MDP(9) + (-g(1) * (-t88 * t105 + t122) - g(2) * (t105 * t90 + t126)) * MDP(11) + (-g(1) * (t88 * t103 + t105 * t125) - g(2) * (-t103 * t90 + t105 * t129)) * MDP(12) + (-g(1) * (-pkin(2) * t88 - qJ(3) * t87 + t124) - g(2) * (pkin(2) * t90 + qJ(3) * t89 + t132)) * MDP(14) + (-g(1) * (pkin(3) * t122 - pkin(4) * t80 - qJ(5) * t79 + t87 * t106 - t88 * t98 + t124) - g(2) * (pkin(3) * t126 + pkin(4) * t84 + qJ(5) * t83 - t89 * t106 + t90 * t98 + t132)) * MDP(25) + (g(1) * t143 - g(2) * t75) * MDP(31) + (g(1) * t142 - g(2) * t74) * MDP(32) + t141 * (-g(1) * t79 + g(2) * t83) - t146 * (-g(1) * t80 + g(2) * t84) - t147 * (g(1) * t87 - g(2) * t89); (-g(1) * (-pkin(2) * t89 + qJ(3) * t90) - g(2) * (-pkin(2) * t87 + qJ(3) * t88) - (pkin(2) * t111 + qJ(3) * t108) * t137) * MDP(14) + (t108 * t137 + t118) * MDP(25) * t106 + (-g(1) * (t110 * t90 - t135 * t89) - g(2) * (t110 * t88 - t135 * t87) - (t108 * t110 + t128 * t99) * t137) * MDP(31) + (-g(1) * (-t107 * t90 - t133 * t89) - g(2) * (-t107 * t88 - t133 * t87) - (-t107 * t108 + t127 * t99) * t137) * MDP(32) + t147 * (g(3) * t130 + t118) + (-t146 * t100 + t141 * t99 - MDP(9) - t105 * MDP(11) + t103 * MDP(12) + (-pkin(4) * t100 - qJ(5) * t99 - t98) * MDP(25)) * t112; (MDP(14) + MDP(25)) * t112; (-g(1) * (-pkin(4) * t83 + qJ(5) * t84) - g(2) * (-pkin(4) * t79 + qJ(5) * t80) - g(3) * (-pkin(4) * t85 + qJ(5) * t86)) * MDP(25) + (-MDP(31) * t107 - MDP(32) * t110 + t141) * (g(1) * t84 + g(2) * t80 + g(3) * t86) + t146 * t115; -t115 * MDP(25); (-g(1) * t74 - g(2) * t142 - g(3) * (t104 * t128 + t85 * t110)) * MDP(31) + (g(1) * t75 + g(2) * t143 - g(3) * (t104 * t127 - t85 * t107)) * MDP(32);];
taug  = t1;
