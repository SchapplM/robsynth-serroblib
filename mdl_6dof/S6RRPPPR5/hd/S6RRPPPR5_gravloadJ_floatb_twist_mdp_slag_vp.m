% Calculate Gravitation load on the joints for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:23:45
% EndTime: 2019-03-09 08:23:47
% DurationCPUTime: 0.70s
% Computational Cost: add. (219->93), mult. (522->129), div. (0->0), fcn. (525->8), ass. (0->48)
t109 = cos(qJ(2));
t106 = sin(qJ(2));
t107 = sin(qJ(1));
t110 = cos(qJ(1));
t120 = g(1) * t110 + g(2) * t107;
t140 = t120 * t106;
t144 = g(3) * t109 - t140;
t143 = MDP(11) - MDP(16) + MDP(21);
t142 = -MDP(12) + MDP(17) + MDP(19);
t141 = MDP(10) - MDP(13) - MDP(15) + MDP(20);
t139 = g(1) * t107;
t136 = g(3) * t106;
t99 = t109 * pkin(2);
t97 = t106 * qJ(3);
t134 = t99 + t97;
t103 = sin(pkin(9));
t133 = qJ(4) * t103;
t104 = cos(pkin(9));
t132 = t104 * t109;
t131 = t106 * t107;
t130 = t106 * t110;
t129 = t107 * t109;
t128 = t109 * t110;
t127 = t110 * t103;
t126 = MDP(18) + MDP(22);
t125 = -pkin(1) - t99;
t124 = -pkin(2) - t133;
t123 = pkin(3) * t132 + t109 * t133 + t134;
t122 = t110 * pkin(1) + pkin(2) * t128 + t107 * pkin(7) + qJ(3) * t130;
t100 = t110 * pkin(7);
t83 = t103 * t129 + t104 * t110;
t84 = t104 * t129 - t127;
t121 = -t84 * pkin(3) - qJ(4) * t83 + t100;
t105 = sin(qJ(6));
t108 = cos(qJ(6));
t118 = t105 * t84 + t108 * t83;
t117 = t105 * t83 - t108 * t84;
t116 = t103 * t108 + t104 * t105;
t115 = t103 * t105 - t104 * t108;
t85 = -t107 * t104 + t109 * t127;
t86 = t107 * t103 + t104 * t128;
t113 = t86 * pkin(3) + t85 * qJ(4) + t122;
t111 = (t125 - t97) * t139;
t92 = qJ(3) * t128;
t89 = qJ(3) * t129;
t74 = t105 * t86 + t108 * t85;
t73 = -t105 * t85 + t108 * t86;
t1 = [t120 * MDP(3) + (-g(1) * t100 - g(2) * t122 - t111) * MDP(14) + (-g(1) * t121 - g(2) * t113 - t111) * MDP(18) + (-g(1) * (-qJ(5) * t84 + t121) - g(2) * (pkin(4) * t130 + t86 * qJ(5) + t113) - ((-pkin(4) - qJ(3)) * t106 + t125) * t139) * MDP(22) + (g(1) * t118 - g(2) * t74) * MDP(28) + (-g(1) * t117 - g(2) * t73) * MDP(29) + (t109 * MDP(9) + MDP(2)) * (-g(2) * t110 + t139) + t143 * (g(1) * t84 - g(2) * t86) + t142 * (g(1) * t83 - g(2) * t85) - t141 * (g(1) * t131 - g(2) * t130); (-g(1) * (-pkin(2) * t130 + t92) - g(2) * (-pkin(2) * t131 + t89) - g(3) * t134) * MDP(14) + (-g(1) * t92 - g(2) * t89 - g(3) * t123 + (pkin(3) * t104 - t124) * t140) * MDP(18) + (-g(1) * (pkin(4) * t128 + t92) - g(2) * (pkin(4) * t129 + t89) - g(3) * (qJ(5) * t132 + t123) + (-g(3) * pkin(4) + t120 * (-(-pkin(3) - qJ(5)) * t104 - t124)) * t106) * MDP(22) + t141 * (t120 * t109 + t136) + (-t116 * MDP(28) + t115 * MDP(29) - t142 * t103 - t143 * t104 - MDP(9)) * t144; -(-MDP(14) - t126) * t144; t126 * (-g(1) * t85 - g(2) * t83 - t103 * t136); (-g(1) * t86 - g(2) * t84 - t104 * t136) * MDP(22); (-g(1) * t73 + g(2) * t117) * MDP(28) + (g(1) * t74 + g(2) * t118) * MDP(29) + (t115 * MDP(28) + t116 * MDP(29)) * t136;];
taug  = t1;
