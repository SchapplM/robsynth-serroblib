% Calculate Gravitation load on the joints for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRPP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:44:47
% EndTime: 2019-03-09 04:44:49
% DurationCPUTime: 0.59s
% Computational Cost: add. (362->90), mult. (490->123), div. (0->0), fcn. (476->8), ass. (0->49)
t152 = MDP(14) - MDP(23) + MDP(28);
t149 = MDP(20) + MDP(22) + MDP(26);
t148 = MDP(21) - MDP(24) - MDP(27);
t112 = pkin(9) + qJ(3);
t109 = sin(t112);
t117 = sin(qJ(1));
t119 = cos(qJ(1));
t97 = g(1) * t119 + g(2) * t117;
t151 = t109 * t97;
t110 = cos(t112);
t150 = -t110 * pkin(3) - t109 * pkin(8);
t147 = -pkin(4) - pkin(5);
t146 = g(1) * t117;
t115 = -pkin(7) - qJ(2);
t144 = g(2) * t115;
t116 = sin(qJ(4));
t142 = qJ(5) * t116;
t141 = t109 * t116;
t118 = cos(qJ(4));
t140 = t109 * t118;
t139 = t109 * t119;
t138 = t110 * t117;
t137 = t110 * t118;
t136 = t110 * t119;
t135 = t117 * t116;
t134 = t117 * t118;
t133 = t118 * t119;
t132 = t119 * t116;
t131 = MDP(25) + MDP(29);
t130 = -pkin(3) - t142;
t89 = t110 * t135 + t133;
t90 = t110 * t134 - t132;
t129 = -t89 * pkin(4) + qJ(5) * t90;
t91 = t110 * t132 - t134;
t92 = t110 * t133 + t135;
t128 = -t91 * pkin(4) + qJ(5) * t92;
t127 = pkin(4) * t137 + t110 * t142 - t150;
t96 = -g(2) * t119 + t146;
t114 = cos(pkin(9));
t106 = pkin(2) * t114 + pkin(1);
t125 = -t106 + t150;
t124 = -t90 * pkin(4) - t89 * qJ(5) - t115 * t119;
t123 = pkin(3) * t136 + t92 * pkin(4) + pkin(8) * t139 + t91 * qJ(5) + t119 * t106;
t121 = g(1) * t91 + g(2) * t89 + g(3) * t141;
t82 = -g(3) * t110 + t151;
t101 = pkin(8) * t136;
t98 = pkin(8) * t138;
t93 = qJ(5) * t140;
t1 = [(-g(1) * (-t117 * pkin(1) + qJ(2) * t119) - g(2) * (pkin(1) * t119 + t117 * qJ(2))) * MDP(7) + (-g(1) * t124 - g(2) * t123 + (-g(1) * t125 + t144) * t117) * MDP(25) + (-g(1) * (-t90 * pkin(5) + t124) - g(2) * (t92 * pkin(5) - qJ(6) * t139 + t123) + (-g(1) * (qJ(6) * t109 + t125) + t144) * t117) * MDP(29) + (MDP(3) - MDP(6)) * t97 - t152 * (-g(2) * t139 + t109 * t146) + t149 * (g(1) * t90 - g(2) * t92) - t148 * (g(1) * t89 - g(2) * t91) + (t110 * MDP(13) + MDP(4) * t114 - MDP(5) * sin(pkin(9)) + MDP(2)) * t96; (-MDP(7) - t131) * t96; (-g(1) * t101 - g(2) * t98 - g(3) * t127 + (pkin(4) * t118 - t130) * t151) * MDP(25) + (-g(1) * (-qJ(6) * t136 + t101) - g(2) * (-qJ(6) * t138 + t98) - g(3) * (pkin(5) * t137 + t127) + (g(3) * qJ(6) + t97 * (-t118 * t147 - t130)) * t109) * MDP(29) + t152 * (g(3) * t109 + t110 * t97) + (-t148 * t116 + t149 * t118 + MDP(13)) * t82; (-g(1) * t128 - g(2) * t129 - g(3) * (-pkin(4) * t141 + t93)) * MDP(25) + (-g(1) * (-pkin(5) * t91 + t128) - g(2) * (-pkin(5) * t89 + t129) - g(3) * (t141 * t147 + t93)) * MDP(29) + t149 * t121 + t148 * (g(1) * t92 + g(2) * t90 + g(3) * t140); -t131 * t121; t82 * MDP(29);];
taug  = t1;
