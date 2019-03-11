% Calculate Gravitation load on the joints for
% S6RRRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:01:20
% EndTime: 2019-03-09 17:01:22
% DurationCPUTime: 0.85s
% Computational Cost: add. (431->135), mult. (792->212), div. (0->0), fcn. (905->12), ass. (0->71)
t144 = sin(qJ(2));
t145 = sin(qJ(1));
t148 = cos(qJ(2));
t149 = cos(qJ(1));
t184 = cos(pkin(6));
t161 = t149 * t184;
t116 = t144 * t145 - t148 * t161;
t142 = sin(qJ(5));
t146 = cos(qJ(5));
t117 = t144 * t161 + t145 * t148;
t138 = qJ(3) + pkin(11);
t135 = sin(t138);
t136 = cos(t138);
t139 = sin(pkin(6));
t171 = t139 * t149;
t99 = t117 * t136 - t135 * t171;
t190 = -t116 * t146 + t142 * t99;
t183 = t116 * t142;
t189 = t146 * t99 + t183;
t187 = MDP(10) - MDP(18);
t185 = g(3) * t139;
t181 = t117 * t142;
t162 = t145 * t184;
t118 = t149 * t144 + t148 * t162;
t180 = t118 * t142;
t119 = -t144 * t162 + t148 * t149;
t179 = t119 * t142;
t143 = sin(qJ(3));
t178 = t119 * t143;
t177 = t136 * t142;
t176 = t136 * t146;
t175 = t139 * t144;
t174 = t139 * t145;
t147 = cos(qJ(3));
t173 = t139 * t147;
t172 = t139 * t148;
t170 = t142 * t148;
t169 = t146 * t148;
t134 = pkin(3) * t147 + pkin(2);
t141 = -qJ(4) - pkin(9);
t168 = -t116 * t134 - t117 * t141;
t167 = -t118 * t134 - t119 * t141;
t166 = t143 * t174;
t165 = t143 * t175;
t164 = t145 * t173;
t163 = t147 * t171;
t127 = t143 * t171;
t160 = t184 * t147;
t159 = t117 * t147 - t127;
t157 = t149 * pkin(1) + pkin(3) * t166 + pkin(8) * t174 - t118 * t141 + t119 * t134;
t103 = t119 * t136 + t135 * t174;
t95 = -t103 * t142 + t118 * t146;
t133 = pkin(5) * t146 + pkin(4);
t140 = -qJ(6) - pkin(10);
t156 = t133 * t136 - t135 * t140;
t98 = t117 * t135 + t136 * t171;
t155 = t117 * t143 + t163;
t154 = -pkin(1) * t145 + pkin(3) * t127 + pkin(8) * t171 + t116 * t141 - t117 * t134;
t102 = t119 * t135 - t136 * t174;
t110 = t135 * t175 - t136 * t184;
t153 = g(1) * t102 + g(2) * t98 + g(3) * t110;
t152 = -g(1) * t118 - g(2) * t116 + g(3) * t172;
t151 = g(1) * t119 + g(2) * t117 + g(3) * t175;
t132 = pkin(3) * t160;
t124 = pkin(3) * t164;
t120 = t134 * t172;
t111 = t135 * t184 + t136 * t175;
t105 = t119 * t147 + t166;
t104 = t164 - t178;
t96 = t103 * t146 + t180;
t1 = [(g(1) * t145 - g(2) * t149) * MDP(2) + (g(1) * t149 + g(2) * t145) * MDP(3) + (g(1) * t117 - g(2) * t119) * MDP(9) + (g(1) * t159 - g(2) * t105) * MDP(16) + (-g(1) * t155 - g(2) * t104) * MDP(17) + (-g(1) * t154 - g(2) * t157) * MDP(19) + (g(1) * t189 - g(2) * t96) * MDP(25) + (-g(1) * t190 - g(2) * t95) * MDP(26) + (g(1) * t98 - g(2) * t102) * MDP(27) + (-g(1) * (-pkin(5) * t183 - t133 * t99 + t140 * t98 + t154) - g(2) * (pkin(5) * t180 - t102 * t140 + t103 * t133 + t157)) * MDP(28) - t187 * (g(1) * t116 - g(2) * t118); (-g(1) * t167 - g(2) * t168 - g(3) * (-t141 * t175 + t120)) * MDP(19) + (-g(1) * (-t118 * t176 + t179) - g(2) * (-t116 * t176 + t181) - (t136 * t169 + t142 * t144) * t185) * MDP(25) + (-g(1) * (t118 * t177 + t119 * t146) - g(2) * (t116 * t177 + t117 * t146) - (-t136 * t170 + t144 * t146) * t185) * MDP(26) + (-g(1) * (pkin(5) * t179 - t118 * t156 + t167) - g(2) * (pkin(5) * t181 - t116 * t156 + t168) - g(3) * t120 - (t156 * t148 + (pkin(5) * t142 - t141) * t144) * t185) * MDP(28) + t187 * t151 + (-t147 * MDP(16) + t143 * MDP(17) - t135 * MDP(27) - MDP(9)) * t152; (-g(1) * t104 + g(2) * t155 - g(3) * (t160 - t165)) * MDP(16) + (g(1) * t105 + g(2) * t159 - g(3) * (-t143 * t184 - t144 * t173)) * MDP(17) + (-g(1) * t124 - g(3) * t132 + (g(2) * t163 + t143 * t151) * pkin(3)) * MDP(19) + (-g(1) * t103 - g(2) * t99 - g(3) * t111) * MDP(27) + (-g(1) * (-pkin(3) * t178 - t102 * t133 - t103 * t140 + t124) - g(2) * (-pkin(3) * t155 - t98 * t133 - t99 * t140) - g(3) * (-pkin(3) * t165 - t110 * t133 - t111 * t140 + t132)) * MDP(28) + (MDP(25) * t146 - MDP(26) * t142) * t153; (MDP(19) + MDP(28)) * t152; (g(1) * t96 + g(2) * t189 - g(3) * (-t111 * t146 + t139 * t170)) * MDP(26) + (pkin(5) * MDP(28) + MDP(25)) * (g(2) * t190 - g(3) * (-t111 * t142 - t139 * t169) - g(1) * t95); -t153 * MDP(28);];
taug  = t1;
