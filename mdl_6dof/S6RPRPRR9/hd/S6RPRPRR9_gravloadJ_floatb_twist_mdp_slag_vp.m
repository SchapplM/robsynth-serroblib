% Calculate Gravitation load on the joints for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRPRR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:05:08
% EndTime: 2019-03-09 04:05:11
% DurationCPUTime: 0.90s
% Computational Cost: add. (591->124), mult. (1601->217), div. (0->0), fcn. (2059->16), ass. (0->71)
t147 = sin(pkin(7));
t145 = sin(pkin(13));
t149 = cos(pkin(13));
t155 = sin(qJ(3));
t159 = cos(qJ(3));
t171 = t159 * t145 + t155 * t149;
t126 = t171 * t147;
t151 = cos(pkin(7));
t128 = t171 * t151;
t152 = cos(pkin(6));
t150 = cos(pkin(12));
t160 = cos(qJ(1));
t177 = t160 * t150;
t146 = sin(pkin(12));
t156 = sin(qJ(1));
t181 = t156 * t146;
t132 = -t152 * t177 + t181;
t178 = t160 * t146;
t180 = t156 * t150;
t133 = t152 * t178 + t180;
t136 = t155 * t145 - t159 * t149;
t148 = sin(pkin(6));
t184 = t148 * t160;
t111 = t126 * t184 + t132 * t128 + t133 * t136;
t120 = -t132 * t147 + t151 * t184;
t154 = sin(qJ(5));
t158 = cos(qJ(5));
t102 = t111 * t158 + t120 * t154;
t153 = sin(qJ(6));
t157 = cos(qJ(6));
t125 = t136 * t147;
t127 = t136 * t151;
t167 = t125 * t184 + t132 * t127 - t133 * t171;
t193 = t102 * t153 - t157 * t167;
t192 = t102 * t157 + t153 * t167;
t191 = t111 * t154 - t120 * t158;
t188 = pkin(9) + qJ(4);
t187 = t147 * t155;
t186 = t148 * t150;
t185 = t148 * t156;
t183 = t151 * t155;
t182 = t153 * t158;
t179 = t157 * t158;
t176 = t160 * pkin(1) + qJ(2) * t185;
t175 = t147 * t184;
t174 = -t156 * pkin(1) + qJ(2) * t184;
t172 = g(1) * t156 - g(2) * t160;
t170 = g(3) * t146 * t148 + g(2) * t133;
t169 = t132 * t151 + t175;
t134 = -t152 * t180 - t178;
t168 = t134 * t151 + t147 * t185;
t164 = t132 * t183 - t133 * t159 + t155 * t175;
t135 = -t152 * t181 + t177;
t163 = t126 * t185 + t134 * t128 - t135 * t136;
t162 = t152 * t126 + (t128 * t150 - t136 * t146) * t148;
t161 = g(2) * t169 - g(3) * (t147 * t152 + t151 * t186);
t143 = t159 * pkin(3) + pkin(2);
t131 = -t147 * t186 + t152 * t151;
t130 = pkin(3) * t183 - t188 * t147;
t129 = pkin(3) * t187 + t188 * t151;
t122 = -t134 * t147 + t151 * t185;
t119 = t135 * t159 + t168 * t155;
t118 = -t135 * t155 + t168 * t159;
t116 = -t152 * t125 + (-t127 * t150 - t146 * t171) * t148;
t113 = -t125 * t185 - t134 * t127 - t135 * t171;
t106 = t131 * t154 + t158 * t162;
t104 = t122 * t154 + t158 * t163;
t103 = t122 * t158 - t154 * t163;
t99 = t104 * t157 - t113 * t153;
t98 = -t104 * t153 - t113 * t157;
t1 = [t172 * MDP(2) + (g(1) * t133 - g(2) * t135) * MDP(4) + (-g(1) * t132 - g(2) * t134) * MDP(5) + (-g(1) * t174 - g(2) * t176) * MDP(7) + (-g(1) * t164 - g(2) * t119) * MDP(13) + (-g(1) * (t133 * t155 + t169 * t159) - g(2) * t118) * MDP(14) + (-g(1) * t120 - g(2) * t122) * MDP(15) + (-g(1) * (t129 * t184 + t132 * t130 - t133 * t143 + t174) - g(2) * (t129 * t185 + t134 * t130 + t135 * t143 + t176)) * MDP(16) + (-g(1) * t102 - g(2) * t104) * MDP(22) + (g(1) * t191 - g(2) * t103) * MDP(23) + (-g(1) * t192 - g(2) * t99) * MDP(29) + (g(1) * t193 - g(2) * t98) * MDP(30) + (-t148 * MDP(6) + MDP(3)) * (g(1) * t160 + g(2) * t156); (MDP(16) + MDP(7)) * (-g(3) * t152 - t172 * t148); (-g(1) * t118 + t170 * t155 + t161 * t159) * MDP(13) + (g(1) * t119 - g(2) * t164 - g(3) * (-t152 * t187 + (-t146 * t159 - t150 * t183) * t148)) * MDP(14) + ((g(1) * t135 + t170) * t155 + (-g(1) * t168 + t161) * t159) * pkin(3) * MDP(16) + (-g(1) * (t113 * t179 + t153 * t163) - g(2) * (-t111 * t153 + t167 * t179) - g(3) * (t116 * t179 + t153 * t162)) * MDP(29) + (-g(1) * (-t113 * t182 + t157 * t163) - g(2) * (-t111 * t157 - t167 * t182) - g(3) * (-t116 * t182 + t157 * t162)) * MDP(30) + (-t158 * MDP(22) + MDP(23) * t154) * (g(1) * t113 + g(2) * t167 + g(3) * t116); (-g(1) * t122 + g(2) * t120 - g(3) * t131) * MDP(16); (g(1) * t104 - g(2) * t102 + g(3) * t106) * MDP(23) + (-MDP(29) * t157 + MDP(30) * t153 - MDP(22)) * (g(1) * t103 + g(2) * t191 + g(3) * (t131 * t158 - t154 * t162)); (-g(1) * t98 - g(2) * t193 - g(3) * (-t106 * t153 - t116 * t157)) * MDP(29) + (g(1) * t99 - g(2) * t192 - g(3) * (-t106 * t157 + t116 * t153)) * MDP(30);];
taug  = t1;
