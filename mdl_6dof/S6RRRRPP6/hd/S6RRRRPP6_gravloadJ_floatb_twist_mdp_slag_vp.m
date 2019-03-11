% Calculate Gravitation load on the joints for
% S6RRRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:15:05
% EndTime: 2019-03-09 21:15:06
% DurationCPUTime: 0.73s
% Computational Cost: add. (513->112), mult. (686->147), div. (0->0), fcn. (693->8), ass. (0->59)
t187 = MDP(24) - MDP(27) - MDP(30);
t183 = MDP(23) - MDP(26) + MDP(31);
t134 = sin(qJ(3));
t136 = sin(qJ(1));
t137 = cos(qJ(3));
t138 = cos(qJ(2));
t139 = cos(qJ(1));
t162 = t138 * t139;
t113 = -t134 * t162 + t136 * t137;
t177 = g(2) * t136;
t152 = g(1) * t139 + t177;
t184 = MDP(10) - MDP(25) - MDP(29);
t181 = pkin(3) * t134;
t133 = qJ(3) + qJ(4);
t129 = cos(t133);
t180 = pkin(4) * t129;
t179 = g(1) * t136;
t135 = sin(qJ(2));
t176 = g(3) * t135;
t174 = -pkin(4) - qJ(6);
t140 = -pkin(9) - pkin(8);
t173 = pkin(5) - t140;
t128 = sin(t133);
t172 = qJ(5) * t128;
t163 = t136 * t138;
t106 = t128 * t163 + t129 * t139;
t171 = qJ(6) * t106;
t161 = t139 * t128;
t108 = -t136 * t129 + t138 * t161;
t170 = qJ(6) * t108;
t169 = t128 * t135;
t168 = t129 * t135;
t167 = t134 * t136;
t166 = t134 * t139;
t165 = t135 * t140;
t127 = pkin(3) * t137 + pkin(2);
t119 = t138 * t127;
t159 = t174 * t128;
t158 = t173 * t139;
t157 = -pkin(1) - t119;
t107 = t129 * t163 - t161;
t156 = -t106 * pkin(4) + qJ(5) * t107;
t109 = t128 * t136 + t129 * t162;
t155 = -t108 * pkin(4) + qJ(5) * t109;
t154 = -t127 - t172;
t153 = g(3) * (t119 + (t172 + t180) * t138);
t144 = g(1) * t109 + g(2) * t107 + g(3) * t168;
t96 = g(1) * t108 + g(2) * t106 + g(3) * t169;
t150 = t187 * t144 + t183 * t96;
t148 = t152 * t138;
t111 = t134 * t163 + t137 * t139;
t147 = pkin(3) * t166 - t107 * pkin(4) + t139 * pkin(7) - qJ(5) * t106 + t136 * t165;
t145 = t139 * pkin(1) + pkin(3) * t167 + t109 * pkin(4) + t136 * pkin(7) + t108 * qJ(5) + t127 * t162;
t143 = t113 * pkin(3) + t155;
t141 = -pkin(3) * t111 + t156;
t117 = qJ(5) * t168;
t114 = t137 * t162 + t167;
t112 = -t137 * t163 + t166;
t1 = [t152 * MDP(3) + (-g(1) * t112 - g(2) * t114) * MDP(16) + (-g(1) * t111 - g(2) * t113) * MDP(17) + (-g(1) * (t136 * t157 + t147) - g(2) * (-t139 * t165 + t145)) * MDP(28) + (-g(1) * (-qJ(6) * t107 + t147) - g(2) * (t109 * qJ(6) + t135 * t158 + t145) - (-t135 * pkin(5) + t157) * t179) * MDP(32) + t183 * (g(1) * t107 - g(2) * t109) - t187 * (g(1) * t106 - g(2) * t108) + (t138 * MDP(9) - t184 * t135 + MDP(2)) * (-g(2) * t139 + t179); (-t153 + t140 * t148 + (g(3) * t140 + t152 * (-t154 + t180)) * t135) * MDP(28) + (-t153 + (-g(3) * qJ(6) * t129 - g(1) * t158 - t173 * t177) * t138 + (-g(3) * t173 + t152 * (-t129 * t174 - t154)) * t135) * MDP(32) + t184 * (t148 + t176) + (t137 * MDP(16) - t134 * MDP(17) - t128 * t187 + t183 * t129 + MDP(9)) * (-g(3) * t138 + t135 * t152); (-g(1) * t113 + g(2) * t111 + t134 * t176) * MDP(16) + (g(1) * t114 - g(2) * t112 + t137 * t176) * MDP(17) + (-g(1) * t143 - g(2) * t141 - g(3) * (t117 + (-pkin(4) * t128 - t181) * t135)) * MDP(28) + (-g(1) * (t143 - t170) - g(2) * (t141 - t171) - g(3) * t117 - (t159 - t181) * t176) * MDP(32) + t150; (-g(1) * t155 - g(2) * t156 - g(3) * (-pkin(4) * t169 + t117)) * MDP(28) + (-g(1) * (t155 - t170) - g(2) * (t156 - t171) - g(3) * (t135 * t159 + t117)) * MDP(32) + t150; -(MDP(28) + MDP(32)) * t96; -t144 * MDP(32);];
taug  = t1;
