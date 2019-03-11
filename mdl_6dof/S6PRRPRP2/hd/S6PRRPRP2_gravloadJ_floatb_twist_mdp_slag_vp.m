% Calculate Gravitation load on the joints for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:32:16
% EndTime: 2019-03-08 21:32:19
% DurationCPUTime: 0.78s
% Computational Cost: add. (457->109), mult. (839->175), div. (0->0), fcn. (991->12), ass. (0->66)
t181 = MDP(19) + MDP(21);
t180 = MDP(20) - MDP(23);
t135 = sin(pkin(10));
t137 = cos(pkin(10));
t144 = cos(qJ(2));
t141 = sin(qJ(2));
t175 = cos(pkin(6));
t155 = t141 * t175;
t121 = t135 * t144 + t137 * t155;
t140 = sin(qJ(3));
t136 = sin(pkin(6));
t143 = cos(qJ(3));
t166 = t136 * t143;
t158 = t137 * t166;
t179 = -t121 * t140 - t158;
t123 = -t135 * t155 + t137 * t144;
t173 = t123 * t140;
t134 = qJ(3) + pkin(11);
t133 = cos(t134);
t139 = sin(qJ(5));
t172 = t133 * t139;
t142 = cos(qJ(5));
t171 = t133 * t142;
t170 = t135 * t136;
t169 = t136 * t137;
t168 = t136 * t140;
t167 = t136 * t141;
t165 = t136 * t144;
t138 = -qJ(4) - pkin(8);
t164 = t138 * t141;
t163 = t142 * t144;
t154 = t144 * t175;
t120 = t135 * t141 - t137 * t154;
t131 = pkin(3) * t143 + pkin(2);
t162 = -t120 * t131 - t121 * t138;
t122 = t135 * t154 + t137 * t141;
t161 = -t122 * t131 - t123 * t138;
t160 = MDP(13) + MDP(24);
t159 = t135 * t166;
t157 = t140 * t167;
t156 = t139 * t165;
t153 = t175 * t143;
t132 = sin(t134);
t152 = pkin(4) * t133 + pkin(9) * t132;
t113 = t132 * t175 + t133 * t167;
t106 = t113 * t139 + t136 * t163;
t103 = t121 * t133 - t132 * t169;
t93 = t103 * t139 - t120 * t142;
t105 = t123 * t133 + t132 * t170;
t95 = t105 * t139 - t122 * t142;
t150 = g(1) * t95 + g(2) * t93 + g(3) * t106;
t146 = -g(1) * t122 - g(2) * t120 + g(3) * t165;
t145 = g(1) * t123 + g(2) * t121 + g(3) * t167;
t130 = pkin(3) * t153;
t125 = pkin(3) * t159;
t124 = t131 * t165;
t109 = (t133 * t163 + t139 * t141) * t136;
t108 = t133 * t156 - t142 * t167;
t107 = t113 * t142 - t156;
t100 = -t122 * t171 + t123 * t139;
t99 = -t122 * t172 - t123 * t142;
t98 = -t120 * t171 + t121 * t139;
t97 = -t120 * t172 - t121 * t142;
t96 = t105 * t142 + t122 * t139;
t94 = t103 * t142 + t120 * t139;
t1 = [(-MDP(1) - t160) * g(3); (-g(1) * t161 - g(2) * t162 - g(3) * (-t136 * t164 + t124)) * MDP(13) + (-g(1) * (pkin(5) * t100 + qJ(6) * t99 - t122 * t152 + t161) - g(2) * (pkin(5) * t98 + qJ(6) * t97 - t120 * t152 + t162) + (-t109 * pkin(5) - t108 * qJ(6) - t124 - (t144 * t152 - t164) * t136) * g(3)) * MDP(24) + t181 * (-g(1) * t100 - g(2) * t98 - g(3) * t109) + t180 * (g(1) * t99 + g(2) * t97 + g(3) * t108) + (MDP(4) - MDP(12)) * t145 + (-MDP(10) * t143 + MDP(11) * t140 - t132 * MDP(22) - MDP(3)) * t146; (-g(1) * (t159 - t173) - g(2) * t179 - g(3) * (t153 - t157)) * MDP(10) + (-g(1) * (-t123 * t143 - t135 * t168) - g(2) * (-t121 * t143 + t137 * t168) - g(3) * (-t140 * t175 - t141 * t166)) * MDP(11) + (-g(1) * t125 - g(3) * t130 + (g(2) * t158 + t140 * t145) * pkin(3)) * MDP(13) + (-g(1) * t105 - g(2) * t103 - g(3) * t113) * MDP(22) + (-g(1) * (-pkin(3) * t173 + pkin(9) * t105 + t125) - g(2) * (t179 * pkin(3) + pkin(9) * t103) - g(3) * (-pkin(3) * t157 + pkin(9) * t113 + t130)) * MDP(24) + ((-pkin(5) * t142 - qJ(6) * t139 - pkin(4)) * MDP(24) - t181 * t142 + t180 * t139) * (g(3) * (-t132 * t167 + t133 * t175) + g(2) * (-t121 * t132 - t133 * t169) + g(1) * (-t123 * t132 + t133 * t170)); t160 * t146; (-g(1) * (-pkin(5) * t95 + qJ(6) * t96) - g(2) * (-pkin(5) * t93 + qJ(6) * t94) - g(3) * (-pkin(5) * t106 + qJ(6) * t107)) * MDP(24) + t181 * t150 + t180 * (g(1) * t96 + g(2) * t94 + g(3) * t107); -t150 * MDP(24);];
taug  = t1;
