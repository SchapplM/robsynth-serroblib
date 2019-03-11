% Calculate Gravitation load on the joints for
% S6RRPRPP3
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
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPRPP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:57:10
% EndTime: 2019-03-09 09:57:12
% DurationCPUTime: 0.78s
% Computational Cost: add. (411->102), mult. (574->137), div. (0->0), fcn. (559->8), ass. (0->52)
t161 = MDP(21) - MDP(24) - MDP(27);
t158 = MDP(20) - MDP(23) + MDP(28);
t120 = cos(qJ(1));
t118 = sin(qJ(1));
t152 = g(2) * t118;
t131 = g(1) * t120 + t152;
t117 = sin(qJ(2));
t119 = cos(qJ(2));
t93 = -g(3) * t119 + t131 * t117;
t156 = MDP(10) - MDP(13) - MDP(22) - MDP(26);
t113 = pkin(9) + qJ(4);
t108 = cos(t113);
t155 = pkin(4) * t108;
t154 = g(1) * t118;
t150 = -pkin(4) - qJ(6);
t116 = -pkin(8) - qJ(3);
t149 = pkin(5) - t116;
t107 = sin(t113);
t148 = qJ(5) * t107;
t147 = t107 * t117;
t146 = t108 * t117;
t114 = sin(pkin(9));
t145 = t114 * t120;
t144 = t116 * t117;
t143 = t118 * t114;
t142 = t118 * t119;
t115 = cos(pkin(9));
t106 = pkin(3) * t115 + pkin(2);
t99 = t119 * t106;
t141 = t119 * t120;
t140 = t120 * t107;
t139 = t120 * pkin(1) + t118 * pkin(7);
t138 = MDP(25) + MDP(29);
t137 = t149 * t120;
t136 = -pkin(1) - t99;
t89 = t107 * t142 + t108 * t120;
t90 = t108 * t142 - t140;
t135 = -t89 * pkin(4) + qJ(5) * t90;
t91 = -t118 * t108 + t119 * t140;
t92 = t118 * t107 + t108 * t141;
t134 = -t91 * pkin(4) + qJ(5) * t92;
t133 = -t106 - t148;
t132 = g(3) * (t99 + (t148 + t155) * t119);
t129 = pkin(2) * t119 + qJ(3) * t117;
t126 = t131 * t119;
t110 = t120 * pkin(7);
t125 = pkin(3) * t145 - t90 * pkin(4) - qJ(5) * t89 + t118 * t144 + t110;
t124 = pkin(3) * t143 + t92 * pkin(4) + t91 * qJ(5) + t106 * t141 + t139;
t122 = g(1) * t91 + g(2) * t89 + g(3) * t147;
t121 = g(1) * t92 + g(2) * t90 + g(3) * t146;
t97 = qJ(5) * t146;
t1 = [t131 * MDP(3) + (-g(1) * (-t115 * t142 + t145) - g(2) * (t115 * t141 + t143)) * MDP(11) + (-g(1) * (t114 * t142 + t115 * t120) - g(2) * (-t114 * t141 + t118 * t115)) * MDP(12) + (-g(1) * t110 - g(2) * (t129 * t120 + t139) - (-pkin(1) - t129) * t154) * MDP(14) + (-g(1) * (t136 * t118 + t125) - g(2) * (-t120 * t144 + t124)) * MDP(25) + (-g(1) * (-qJ(6) * t90 + t125) - g(2) * (t92 * qJ(6) + t117 * t137 + t124) - (-pkin(5) * t117 + t136) * t154) * MDP(29) + t158 * (g(1) * t90 - g(2) * t92) - t161 * (g(1) * t89 - g(2) * t91) + (t119 * MDP(9) - t156 * t117 + MDP(2)) * (-g(2) * t120 + t154); (-g(3) * t129 + t131 * (pkin(2) * t117 - qJ(3) * t119)) * MDP(14) + (-t132 + t116 * t126 + (g(3) * t116 + t131 * (-t133 + t155)) * t117) * MDP(25) + (-t132 + (-g(3) * qJ(6) * t108 - g(1) * t137 - t149 * t152) * t119 + (-g(3) * t149 + t131 * (-t150 * t108 - t133)) * t117) * MDP(29) + t156 * (g(3) * t117 + t126) + (t115 * MDP(11) - t114 * MDP(12) - t107 * t161 + t108 * t158 + MDP(9)) * t93; (-MDP(14) - t138) * t93; (-g(1) * t134 - g(2) * t135 - g(3) * (-pkin(4) * t147 + t97)) * MDP(25) + (-g(1) * (-qJ(6) * t91 + t134) - g(2) * (-qJ(6) * t89 + t135) - g(3) * (t150 * t147 + t97)) * MDP(29) + t158 * t122 + t161 * t121; -t138 * t122; -t121 * MDP(29);];
taug  = t1;
