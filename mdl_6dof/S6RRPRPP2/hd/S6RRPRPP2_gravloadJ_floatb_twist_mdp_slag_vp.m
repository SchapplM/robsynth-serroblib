% Calculate Gravitation load on the joints for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPRPP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:52:29
% EndTime: 2019-03-09 09:52:31
% DurationCPUTime: 0.62s
% Computational Cost: add. (361->95), mult. (501->128), div. (0->0), fcn. (485->8), ass. (0->53)
t159 = MDP(21) - MDP(26);
t156 = MDP(18) + MDP(20) + MDP(24);
t155 = MDP(19) - MDP(22) - MDP(25);
t114 = qJ(2) + pkin(9);
t110 = sin(t114);
t118 = sin(qJ(1));
t121 = cos(qJ(1));
t98 = g(1) * t121 + g(2) * t118;
t158 = t110 * t98;
t111 = cos(t114);
t157 = -t111 * pkin(3) - t110 * pkin(8);
t154 = -pkin(4) - pkin(5);
t117 = sin(qJ(2));
t153 = pkin(2) * t117;
t115 = -qJ(3) - pkin(7);
t151 = g(2) * t115;
t116 = sin(qJ(4));
t149 = qJ(5) * t116;
t148 = t110 * t116;
t119 = cos(qJ(4));
t147 = t110 * t119;
t146 = t110 * t121;
t145 = t111 * t118;
t144 = t111 * t119;
t143 = t111 * t121;
t142 = t115 * t121;
t141 = t118 * t116;
t140 = t118 * t119;
t139 = t119 * t121;
t138 = t121 * t116;
t137 = MDP(23) + MDP(27);
t136 = -pkin(3) - t149;
t91 = t111 * t141 + t139;
t92 = t111 * t140 - t138;
t135 = -t91 * pkin(4) + qJ(5) * t92;
t93 = t111 * t138 - t140;
t94 = t111 * t139 + t141;
t134 = -t93 * pkin(4) + qJ(5) * t94;
t133 = pkin(8) * t145 - t118 * t153;
t132 = pkin(8) * t143 - t121 * t153;
t97 = g(1) * t118 - g(2) * t121;
t120 = cos(qJ(2));
t112 = t120 * pkin(2);
t131 = pkin(4) * t144 + t111 * t149 + t112 - t157;
t109 = t112 + pkin(1);
t129 = -t109 + t157;
t128 = -t92 * pkin(4) - t91 * qJ(5) - t142;
t104 = t121 * t109;
t127 = pkin(3) * t143 + t94 * pkin(4) + pkin(8) * t146 + t93 * qJ(5) + t104;
t125 = g(1) * t93 + g(2) * t91 + g(3) * t148;
t123 = -g(3) * t111 + t158;
t95 = qJ(5) * t147;
t1 = [(-g(1) * (-t118 * t109 - t142) - g(2) * (-t118 * t115 + t104)) * MDP(12) + (-g(1) * t128 - g(2) * t127 + (-g(1) * t129 + t151) * t118) * MDP(23) + (-g(1) * (-t92 * pkin(5) + t128) - g(2) * (t94 * pkin(5) - qJ(6) * t146 + t127) + (-g(1) * (qJ(6) * t110 + t129) + t151) * t118) * MDP(27) + (MDP(3) - MDP(11)) * t98 + t156 * (g(1) * t92 - g(2) * t94) - t155 * (g(1) * t91 - g(2) * t93) + (-t117 * MDP(10) + t120 * MDP(9) + t159 * t110 + MDP(2)) * t97; (g(3) * t117 + t120 * t98) * MDP(10) + (-g(1) * t132 - g(2) * t133 - g(3) * t131 + (pkin(4) * t119 - t136) * t158) * MDP(23) + (-g(1) * (-qJ(6) * t143 + t132) - g(2) * (-qJ(6) * t145 + t133) - g(3) * (pkin(5) * t144 + t131) + (g(3) * qJ(6) + t98 * (-t119 * t154 - t136)) * t110) * MDP(27) - t159 * (g(3) * t110 + t111 * t98) + (pkin(2) * MDP(12) + MDP(9)) * (-g(3) * t120 + t117 * t98) + (-t155 * t116 + t156 * t119) * t123; (-MDP(12) - t137) * t97; (-g(1) * t134 - g(2) * t135 - g(3) * (-pkin(4) * t148 + t95)) * MDP(23) + (-g(1) * (-pkin(5) * t93 + t134) - g(2) * (-pkin(5) * t91 + t135) - g(3) * (t154 * t148 + t95)) * MDP(27) + t156 * t125 + t155 * (g(1) * t94 + g(2) * t92 + g(3) * t147); -t137 * t125; t123 * MDP(27);];
taug  = t1;
