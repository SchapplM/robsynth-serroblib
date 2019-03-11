% Calculate Gravitation load on the joints for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:53:22
% EndTime: 2019-03-08 23:53:24
% DurationCPUTime: 1.04s
% Computational Cost: add. (658->134), mult. (1835->229), div. (0->0), fcn. (2331->14), ass. (0->67)
t190 = MDP(17) - MDP(20);
t136 = sin(pkin(7));
t142 = sin(qJ(2));
t145 = cos(qJ(2));
t138 = cos(pkin(12));
t179 = cos(pkin(6));
t165 = t138 * t179;
t177 = sin(pkin(12));
t151 = t177 * t142 - t145 * t165;
t137 = sin(pkin(6));
t174 = t137 * t138;
t178 = cos(pkin(7));
t189 = t136 * t174 + t151 * t178;
t161 = t179 * t177;
t150 = t138 * t142 + t145 * t161;
t166 = t137 * t177;
t188 = -t136 * t166 + t150 * t178;
t187 = MDP(11) - MDP(19);
t185 = MDP(18) - MDP(21);
t184 = cos(qJ(3));
t180 = t136 * pkin(9);
t140 = sin(qJ(4));
t176 = t136 * t140;
t144 = cos(qJ(4));
t175 = t136 * t144;
t173 = t137 * t142;
t172 = t137 * t145;
t139 = sin(qJ(6));
t171 = t139 * t140;
t143 = cos(qJ(6));
t170 = t140 * t143;
t168 = t136 * t173;
t167 = t136 * t179;
t141 = sin(qJ(3));
t164 = t141 * t178;
t162 = t178 * t184;
t120 = t141 * t167 + (t184 * t142 + t145 * t164) * t137;
t152 = -t136 * t172 + t179 * t178;
t109 = t120 * t140 - t152 * t144;
t128 = t142 * t165 + t177 * t145;
t106 = t128 * t184 - t189 * t141;
t147 = t151 * t136 - t178 * t174;
t97 = t106 * t140 - t147 * t144;
t129 = t138 * t145 - t142 * t161;
t108 = t129 * t184 - t188 * t141;
t146 = t150 * t136 + t178 * t166;
t99 = t108 * t140 - t146 * t144;
t159 = g(1) * t99 + g(2) * t97 + g(3) * t109;
t126 = (-t142 * t164 + t184 * t145) * t137;
t125 = (t141 * t145 + t142 * t162) * t137;
t119 = t141 * t173 - t162 * t172 - t184 * t167;
t116 = t126 * t144 + t140 * t168;
t115 = t126 * t140 - t144 * t168;
t114 = -t129 * t164 - t150 * t184;
t113 = t129 * t162 - t150 * t141;
t112 = -t128 * t164 - t151 * t184;
t111 = t128 * t162 - t151 * t141;
t110 = t120 * t144 + t152 * t140;
t107 = t129 * t141 + t188 * t184;
t105 = t128 * t141 + t189 * t184;
t104 = t114 * t144 + t129 * t176;
t103 = t114 * t140 - t129 * t175;
t102 = t112 * t144 + t128 * t176;
t101 = t112 * t140 - t128 * t175;
t100 = t108 * t144 + t146 * t140;
t98 = t106 * t144 + t147 * t140;
t1 = [(-MDP(1) - MDP(22)) * g(3); (g(1) * t150 + g(2) * t151 - g(3) * t172) * MDP(3) + (g(1) * t129 + g(2) * t128 + g(3) * t173) * MDP(4) + (-g(1) * t114 - g(2) * t112 - g(3) * t126) * MDP(10) + (-g(1) * (-pkin(2) * t150 + t114 * pkin(3) + t104 * pkin(4) + t113 * pkin(10) + t103 * qJ(5) + t129 * t180) - g(2) * (-pkin(2) * t151 + t112 * pkin(3) + t102 * pkin(4) + t111 * pkin(10) + t101 * qJ(5) + t128 * t180) - g(3) * (t126 * pkin(3) + t116 * pkin(4) + t125 * pkin(10) + t115 * qJ(5) + (pkin(2) * t145 + t142 * t180) * t137)) * MDP(22) + (-g(1) * (t103 * t139 + t113 * t143) - g(2) * (t101 * t139 + t111 * t143) - g(3) * (t115 * t139 + t125 * t143)) * MDP(28) + (-g(1) * (t103 * t143 - t113 * t139) - g(2) * (t101 * t143 - t111 * t139) - g(3) * (t115 * t143 - t125 * t139)) * MDP(29) + t185 * (g(1) * t103 + g(2) * t101 + g(3) * t115) - t190 * (g(1) * t104 + g(2) * t102 + g(3) * t116) + t187 * (g(1) * t113 + g(2) * t111 + g(3) * t125); (-g(1) * (-t107 * t171 + t108 * t143) - g(2) * (-t105 * t171 + t106 * t143) - g(3) * (-t119 * t171 + t120 * t143)) * MDP(28) + (-g(1) * (-t107 * t170 - t108 * t139) - g(2) * (-t105 * t170 - t106 * t139) - g(3) * (-t119 * t170 - t120 * t139)) * MDP(29) + (-pkin(10) * MDP(22) + t187) * (g(1) * t108 + g(2) * t106 + g(3) * t120) + (MDP(10) + t190 * t144 - t185 * t140 + MDP(22) * (pkin(4) * t144 + qJ(5) * t140 + pkin(3))) * (g(1) * t107 + g(2) * t105 + g(3) * t119); (-g(1) * (-pkin(4) * t99 + qJ(5) * t100) - g(2) * (-pkin(4) * t97 + qJ(5) * t98) - g(3) * (-pkin(4) * t109 + qJ(5) * t110)) * MDP(22) + (-MDP(28) * t139 - MDP(29) * t143 + t185) * (g(1) * t100 + g(2) * t98 + g(3) * t110) + t190 * t159; -t159 * MDP(22); (-g(1) * (-t107 * t139 + t143 * t99) - g(2) * (-t105 * t139 + t143 * t97) - g(3) * (t109 * t143 - t119 * t139)) * MDP(28) + (-g(1) * (-t107 * t143 - t139 * t99) - g(2) * (-t105 * t143 - t139 * t97) - g(3) * (-t109 * t139 - t119 * t143)) * MDP(29);];
taug  = t1;
