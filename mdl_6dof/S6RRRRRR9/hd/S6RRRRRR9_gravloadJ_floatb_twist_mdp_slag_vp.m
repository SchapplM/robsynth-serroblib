% Calculate Gravitation load on the joints for
% S6RRRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:39:16
% EndTime: 2019-03-10 05:39:20
% DurationCPUTime: 1.27s
% Computational Cost: add. (870->164), mult. (2250->290), div. (0->0), fcn. (2908->16), ass. (0->72)
t146 = sin(qJ(2));
t149 = cos(qJ(2));
t150 = cos(qJ(1));
t174 = cos(pkin(6));
t159 = t150 * t174;
t176 = sin(qJ(1));
t130 = t176 * t146 - t149 * t159;
t131 = t146 * t159 + t176 * t149;
t145 = sin(qJ(3));
t173 = cos(pkin(7));
t160 = t145 * t173;
t141 = sin(pkin(7));
t142 = sin(pkin(6));
t170 = t141 * t142;
t162 = t150 * t170;
t177 = cos(qJ(3));
t110 = -t130 * t160 + t131 * t177 - t145 * t162;
t161 = t142 * t173;
t123 = t130 * t141 - t150 * t161;
t144 = sin(qJ(4));
t148 = cos(qJ(4));
t101 = t110 * t148 + t123 * t144;
t155 = t173 * t177;
t109 = t130 * t155 + t131 * t145 + t177 * t162;
t140 = qJ(5) + qJ(6);
t138 = sin(t140);
t139 = cos(t140);
t187 = t101 * t138 - t109 * t139;
t186 = t101 * t139 + t109 * t138;
t143 = sin(qJ(5));
t147 = cos(qJ(5));
t185 = t101 * t143 - t109 * t147;
t184 = t101 * t147 + t109 * t143;
t179 = t110 * t144 - t123 * t148;
t156 = t174 * t176;
t152 = t150 * t146 + t149 * t156;
t178 = t152 * t173 - t176 * t170;
t158 = t174 * t141;
t121 = t145 * t158 + (t177 * t146 + t149 * t160) * t142;
t166 = t142 * t149;
t129 = -t141 * t166 + t174 * t173;
t108 = t121 * t148 + t129 * t144;
t167 = t142 * t146;
t120 = t145 * t167 - t155 * t166 - t177 * t158;
t132 = -t146 * t156 + t150 * t149;
t114 = t132 * t177 - t178 * t145;
t124 = t152 * t141 + t176 * t161;
t104 = t114 * t148 + t124 * t144;
t113 = t132 * t145 + t178 * t177;
t96 = -t104 * t138 + t113 * t139;
t97 = t104 * t139 + t113 * t138;
t175 = (-g(1) * t96 + g(2) * t187 - g(3) * (-t108 * t138 + t120 * t139)) * MDP(37) + (g(1) * t97 + g(2) * t186 - g(3) * (-t108 * t139 - t120 * t138)) * MDP(38);
t172 = t138 * t148;
t171 = t139 * t148;
t169 = t141 * t144;
t168 = t141 * t148;
t165 = t143 * t148;
t164 = t147 * t148;
t163 = t141 * t167;
t128 = (-t146 * t160 + t177 * t149) * t142;
t127 = (t145 * t149 + t146 * t155) * t142;
t119 = t128 * t148 + t144 * t163;
t118 = -t132 * t160 - t152 * t177;
t117 = t132 * t155 - t152 * t145;
t116 = -t130 * t177 - t131 * t160;
t115 = -t130 * t145 + t131 * t155;
t106 = t118 * t148 + t132 * t169;
t105 = t116 * t148 + t131 * t169;
t103 = -t114 * t144 + t124 * t148;
t99 = t104 * t147 + t113 * t143;
t98 = -t104 * t143 + t113 * t147;
t1 = [(g(1) * t176 - g(2) * t150) * MDP(2) + (g(1) * t150 + g(2) * t176) * MDP(3) + (g(1) * t131 - g(2) * t132) * MDP(9) + (-g(1) * t130 + g(2) * t152) * MDP(10) + (g(1) * t110 - g(2) * t114) * MDP(16) + (-g(1) * t109 + g(2) * t113) * MDP(17) + (g(1) * t101 - g(2) * t104) * MDP(23) + (-g(1) * t179 - g(2) * t103) * MDP(24) + (g(1) * t184 - g(2) * t99) * MDP(30) + (-g(1) * t185 - g(2) * t98) * MDP(31) + (g(1) * t186 - g(2) * t97) * MDP(37) + (-g(1) * t187 - g(2) * t96) * MDP(38); (g(1) * t152 + g(2) * t130 - g(3) * t166) * MDP(9) + (g(1) * t132 + g(2) * t131 + g(3) * t167) * MDP(10) + (-g(1) * t118 - g(2) * t116 - g(3) * t128) * MDP(16) + (g(1) * t117 + g(2) * t115 + g(3) * t127) * MDP(17) + (-g(1) * t106 - g(2) * t105 - g(3) * t119) * MDP(23) + (-g(1) * (-t118 * t144 + t132 * t168) - g(2) * (-t116 * t144 + t131 * t168) - g(3) * (-t128 * t144 + t148 * t163)) * MDP(24) + (-g(1) * (t106 * t147 + t117 * t143) - g(2) * (t105 * t147 + t115 * t143) - g(3) * (t119 * t147 + t127 * t143)) * MDP(30) + (-g(1) * (-t106 * t143 + t117 * t147) - g(2) * (-t105 * t143 + t115 * t147) - g(3) * (-t119 * t143 + t127 * t147)) * MDP(31) + (-g(1) * (t106 * t139 + t117 * t138) - g(2) * (t105 * t139 + t115 * t138) - g(3) * (t119 * t139 + t127 * t138)) * MDP(37) + (-g(1) * (-t106 * t138 + t117 * t139) - g(2) * (-t105 * t138 + t115 * t139) - g(3) * (-t119 * t138 + t127 * t139)) * MDP(38); (g(1) * t114 + g(2) * t110 + g(3) * t121) * MDP(17) + (-g(1) * (-t113 * t164 + t114 * t143) - g(2) * (-t109 * t164 + t110 * t143) - g(3) * (-t120 * t164 + t121 * t143)) * MDP(30) + (-g(1) * (t113 * t165 + t114 * t147) - g(2) * (t109 * t165 + t110 * t147) - g(3) * (t120 * t165 + t121 * t147)) * MDP(31) + (-g(1) * (-t113 * t171 + t114 * t138) - g(2) * (-t109 * t171 + t110 * t138) - g(3) * (-t120 * t171 + t121 * t138)) * MDP(37) + (-g(1) * (t113 * t172 + t114 * t139) - g(2) * (t109 * t172 + t110 * t139) - g(3) * (t120 * t172 + t121 * t139)) * MDP(38) + (t148 * MDP(23) - MDP(24) * t144 + MDP(16)) * (g(1) * t113 + g(2) * t109 + g(3) * t120); (g(1) * t104 + g(2) * t101 + g(3) * t108) * MDP(24) + (-MDP(30) * t147 + MDP(31) * t143 - MDP(37) * t139 + MDP(38) * t138 - MDP(23)) * (g(1) * t103 - g(2) * t179 + g(3) * (-t121 * t144 + t129 * t148)); (-g(1) * t98 + g(2) * t185 - g(3) * (-t108 * t143 + t120 * t147)) * MDP(30) + (g(1) * t99 + g(2) * t184 - g(3) * (-t108 * t147 - t120 * t143)) * MDP(31) + t175; t175;];
taug  = t1;
