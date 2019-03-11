% Calculate Gravitation load on the joints for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR14_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR14_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:20:37
% EndTime: 2019-03-09 20:20:40
% DurationCPUTime: 1.08s
% Computational Cost: add. (409->121), mult. (908->199), div. (0->0), fcn. (1089->12), ass. (0->53)
t151 = MDP(16) - MDP(19);
t150 = MDP(10) - MDP(18);
t143 = MDP(17) - MDP(20);
t102 = qJ(5) + qJ(6);
t100 = sin(t102);
t101 = cos(t102);
t105 = sin(qJ(3));
t109 = cos(qJ(3));
t103 = sin(pkin(6));
t140 = cos(qJ(1));
t122 = t103 * t140;
t106 = sin(qJ(2));
t107 = sin(qJ(1));
t110 = cos(qJ(2));
t133 = cos(pkin(6));
t116 = t133 * t140;
t92 = t106 * t116 + t107 * t110;
t83 = t92 * t105 + t109 * t122;
t91 = t106 * t107 - t110 * t116;
t148 = t100 * t83 + t101 * t91;
t147 = -t100 * t91 + t101 * t83;
t104 = sin(qJ(5));
t108 = cos(qJ(5));
t146 = t104 * t83 + t108 * t91;
t145 = -t104 * t91 + t108 * t83;
t121 = t107 * t133;
t93 = t140 * t106 + t110 * t121;
t144 = -g(1) * t93 - g(2) * t91;
t139 = g(3) * t103;
t127 = t103 * t110;
t128 = t103 * t109;
t94 = -t106 * t121 + t140 * t110;
t87 = t105 * t94 - t107 * t128;
t77 = -t100 * t93 + t101 * t87;
t78 = t100 * t87 + t101 * t93;
t130 = t103 * t106;
t89 = t105 * t130 - t133 * t109;
t138 = (-g(1) * t77 - g(2) * t147 - g(3) * (t100 * t127 + t89 * t101)) * MDP(34) + (g(1) * t78 + g(2) * t148 - g(3) * (-t89 * t100 + t101 * t127)) * MDP(35);
t132 = t100 * t105;
t131 = t101 * t105;
t129 = t103 * t107;
t126 = t104 * t105;
t125 = t105 * t108;
t124 = t105 * t110;
t123 = t108 * t110;
t84 = -t105 * t122 + t109 * t92;
t117 = -g(1) * t94 - g(2) * t92;
t114 = g(1) * t87 + g(2) * t83 + g(3) * t89;
t90 = t133 * t105 + t106 * t128;
t88 = t105 * t129 + t109 * t94;
t80 = t104 * t87 + t108 * t93;
t79 = -t104 * t93 + t108 * t87;
t1 = [(g(1) * t107 - g(2) * t140) * MDP(2) + (g(1) * t140 + g(2) * t107) * MDP(3) + (g(1) * t92 - g(2) * t94) * MDP(9) + (-g(1) * (-t107 * pkin(1) - t92 * pkin(2) - pkin(3) * t84 + pkin(8) * t122 - t91 * pkin(9) - qJ(4) * t83) - g(2) * (t140 * pkin(1) + t94 * pkin(2) + t88 * pkin(3) + pkin(8) * t129 + t93 * pkin(9) + t87 * qJ(4))) * MDP(21) + (g(1) * t146 - g(2) * t80) * MDP(27) + (g(1) * t145 - g(2) * t79) * MDP(28) + (g(1) * t148 - g(2) * t78) * MDP(34) + (g(1) * t147 - g(2) * t77) * MDP(35) + t143 * (-g(1) * t83 + g(2) * t87) - t151 * (-g(1) * t84 + g(2) * t88) - t150 * (g(1) * t91 - g(2) * t93); ((-t106 * t139 + t117) * pkin(9) + (-t110 * t139 - t144) * (pkin(3) * t109 + qJ(4) * t105 + pkin(2))) * MDP(21) + (-g(1) * (t108 * t94 - t93 * t126) - g(2) * (t108 * t92 - t91 * t126) - (t104 * t124 + t106 * t108) * t139) * MDP(27) + (-g(1) * (-t104 * t94 - t93 * t125) - g(2) * (-t104 * t92 - t91 * t125) - (-t104 * t106 + t105 * t123) * t139) * MDP(28) + (-g(1) * (t101 * t94 - t93 * t132) - g(2) * (t101 * t92 - t91 * t132) - (t100 * t124 + t101 * t106) * t139) * MDP(34) + (-g(1) * (-t100 * t94 - t93 * t131) - g(2) * (-t100 * t92 - t91 * t131) - (-t100 * t106 + t101 * t124) * t139) * MDP(35) + t150 * (g(3) * t130 - t117) + (t105 * t143 - t109 * t151 - MDP(9)) * (g(3) * t127 + t144); (-g(1) * (-pkin(3) * t87 + qJ(4) * t88) - g(2) * (-pkin(3) * t83 + qJ(4) * t84) - g(3) * (-pkin(3) * t89 + qJ(4) * t90)) * MDP(21) + (-MDP(27) * t104 - MDP(28) * t108 - MDP(34) * t100 - MDP(35) * t101 + t143) * (g(1) * t88 + g(2) * t84 + g(3) * t90) + t151 * t114; -t114 * MDP(21); (-g(1) * t79 - g(2) * t145 - g(3) * (t104 * t127 + t89 * t108)) * MDP(27) + (g(1) * t80 + g(2) * t146 - g(3) * (t103 * t123 - t89 * t104)) * MDP(28) + t138; t138;];
taug  = t1;
