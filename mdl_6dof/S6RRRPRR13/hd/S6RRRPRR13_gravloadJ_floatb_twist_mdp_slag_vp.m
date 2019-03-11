% Calculate Gravitation load on the joints for
% S6RRRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR13_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR13_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:05:12
% EndTime: 2019-03-09 20:05:17
% DurationCPUTime: 1.52s
% Computational Cost: add. (788->176), mult. (1932->305), div. (0->0), fcn. (2448->16), ass. (0->72)
t147 = sin(qJ(2));
t149 = cos(qJ(2));
t150 = cos(qJ(1));
t179 = cos(pkin(6));
t161 = t150 * t179;
t181 = sin(qJ(1));
t130 = t147 * t181 - t149 * t161;
t131 = t147 * t161 + t149 * t181;
t146 = sin(qJ(3));
t178 = cos(pkin(7));
t182 = cos(qJ(3));
t158 = t178 * t182;
t142 = sin(pkin(7));
t143 = sin(pkin(6));
t169 = t143 * t150;
t166 = t142 * t169;
t106 = t130 * t158 + t131 * t146 + t166 * t182;
t145 = sin(qJ(6));
t148 = cos(qJ(6));
t162 = t146 * t178;
t107 = -t130 * t162 + t131 * t182 - t146 * t166;
t163 = t143 * t178;
t120 = t130 * t142 - t150 * t163;
t140 = pkin(13) + qJ(5);
t138 = sin(t140);
t139 = cos(t140);
t98 = t107 * t139 + t120 * t138;
t189 = -t106 * t148 + t145 * t98;
t188 = t106 * t145 + t148 * t98;
t185 = t107 * t138 - t120 * t139;
t159 = t179 * t181;
t152 = t150 * t147 + t149 * t159;
t165 = t143 * t181;
t184 = -t142 * t165 + t152 * t178;
t183 = MDP(17) - MDP(20);
t180 = pkin(10) * t142;
t177 = t138 * t142;
t176 = t139 * t142;
t175 = t139 * t145;
t174 = t139 * t148;
t141 = sin(pkin(13));
t173 = t141 * t142;
t144 = cos(pkin(13));
t172 = t142 * t144;
t171 = t143 * t147;
t170 = t143 * t149;
t168 = t147 * t142;
t167 = t143 * t168;
t164 = t142 * t179;
t132 = -t147 * t159 + t149 * t150;
t110 = t132 * t146 + t182 * t184;
t117 = t146 * t171 - t158 * t170 - t164 * t182;
t155 = g(1) * t110 + g(2) * t106 + g(3) * t117;
t121 = t142 * t152 + t163 * t181;
t129 = -t142 * t170 + t178 * t179;
t128 = (-t147 * t162 + t149 * t182) * t143;
t127 = (t146 * t149 + t147 * t158) * t143;
t118 = t146 * t164 + (t147 * t182 + t149 * t162) * t143;
t116 = -t132 * t162 - t152 * t182;
t115 = t132 * t158 - t146 * t152;
t114 = -t130 * t182 - t131 * t162;
t113 = -t130 * t146 + t131 * t158;
t112 = t128 * t139 + t138 * t167;
t111 = t132 * t182 - t146 * t184;
t105 = t118 * t139 + t129 * t138;
t103 = t116 * t139 + t132 * t177;
t102 = t114 * t139 + t131 * t177;
t101 = t111 * t139 + t121 * t138;
t100 = -t111 * t138 + t121 * t139;
t96 = t101 * t148 + t110 * t145;
t95 = -t101 * t145 + t110 * t148;
t1 = [(g(1) * t181 - g(2) * t150) * MDP(2) + (g(1) * t150 + g(2) * t181) * MDP(3) + (g(1) * t131 - g(2) * t132) * MDP(9) + (-g(1) * t130 + g(2) * t152) * MDP(10) + (g(1) * t107 - g(2) * t111) * MDP(16) + (-g(1) * (-t107 * t144 - t120 * t141) - g(2) * (t111 * t144 + t121 * t141)) * MDP(18) + (-g(1) * (t107 * t141 - t120 * t144) - g(2) * (-t111 * t141 + t121 * t144)) * MDP(19) + (-g(1) * (-pkin(1) * t181 - t131 * pkin(2) - pkin(3) * t107 + pkin(9) * t169 - qJ(4) * t106) - g(2) * (t150 * pkin(1) + t132 * pkin(2) + t111 * pkin(3) + pkin(9) * t165 + t110 * qJ(4)) + (g(1) * t120 - g(2) * t121) * pkin(10)) * MDP(21) + (g(1) * t98 - g(2) * t101) * MDP(27) + (-g(1) * t185 - g(2) * t100) * MDP(28) + (g(1) * t188 - g(2) * t96) * MDP(34) + (-g(1) * t189 - g(2) * t95) * MDP(35) + t183 * (-g(1) * t106 + g(2) * t110); (g(1) * t152 + g(2) * t130 - g(3) * t170) * MDP(9) + (g(1) * t132 + g(2) * t131 + g(3) * t171) * MDP(10) + (-g(1) * t116 - g(2) * t114 - g(3) * t128) * MDP(16) + (-g(1) * (t116 * t144 + t132 * t173) - g(2) * (t114 * t144 + t131 * t173) - g(3) * (t128 * t144 + t141 * t167)) * MDP(18) + (-g(1) * (-t116 * t141 + t132 * t172) - g(2) * (-t114 * t141 + t131 * t172) - g(3) * (-t128 * t141 + t144 * t167)) * MDP(19) + (-g(1) * (-pkin(2) * t152 + t116 * pkin(3) + t115 * qJ(4) + t132 * t180) - g(2) * (-pkin(2) * t130 + pkin(3) * t114 + qJ(4) * t113 + t131 * t180) - g(3) * (pkin(3) * t128 + qJ(4) * t127 + (pkin(2) * t149 + pkin(10) * t168) * t143)) * MDP(21) + (-g(1) * t103 - g(2) * t102 - g(3) * t112) * MDP(27) + (-g(1) * (-t116 * t138 + t132 * t176) - g(2) * (-t114 * t138 + t131 * t176) - g(3) * (-t128 * t138 + t139 * t167)) * MDP(28) + (-g(1) * (t103 * t148 + t115 * t145) - g(2) * (t102 * t148 + t113 * t145) - g(3) * (t112 * t148 + t127 * t145)) * MDP(34) + (-g(1) * (-t103 * t145 + t115 * t148) - g(2) * (-t102 * t145 + t113 * t148) - g(3) * (-t112 * t145 + t127 * t148)) * MDP(35) + t183 * (g(1) * t115 + g(2) * t113 + g(3) * t127); (-g(1) * (-pkin(3) * t110 + qJ(4) * t111) - g(2) * (-pkin(3) * t106 + qJ(4) * t107) - g(3) * (-pkin(3) * t117 + qJ(4) * t118)) * MDP(21) + (-g(1) * (-t110 * t174 + t111 * t145) - g(2) * (-t106 * t174 + t107 * t145) - g(3) * (-t117 * t174 + t118 * t145)) * MDP(34) + (-g(1) * (t110 * t175 + t111 * t148) - g(2) * (t106 * t175 + t107 * t148) - g(3) * (t117 * t175 + t118 * t148)) * MDP(35) + t183 * (g(1) * t111 + g(2) * t107 + g(3) * t118) + (MDP(18) * t144 - MDP(19) * t141 + MDP(27) * t139 - MDP(28) * t138 + MDP(16)) * t155; -t155 * MDP(21); (g(1) * t101 + g(2) * t98 + g(3) * t105) * MDP(28) + (-MDP(34) * t148 + MDP(35) * t145 - MDP(27)) * (g(1) * t100 - g(2) * t185 + g(3) * (-t118 * t138 + t129 * t139)); (-g(1) * t95 + g(2) * t189 - g(3) * (-t105 * t145 + t117 * t148)) * MDP(34) + (g(1) * t96 + g(2) * t188 - g(3) * (-t105 * t148 - t117 * t145)) * MDP(35);];
taug  = t1;
