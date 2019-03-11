% Calculate Gravitation load on the joints for
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:42:39
% EndTime: 2019-03-09 12:42:42
% DurationCPUTime: 1.15s
% Computational Cost: add. (652->134), mult. (1157->204), div. (0->0), fcn. (1384->12), ass. (0->59)
t200 = MDP(10) - MDP(13);
t198 = MDP(21) - MDP(30);
t179 = MDP(27) + MDP(29);
t197 = MDP(28) - MDP(31);
t155 = sin(qJ(2));
t156 = sin(qJ(1));
t158 = cos(qJ(2));
t189 = cos(pkin(6));
t196 = cos(qJ(1));
t171 = t189 * t196;
t133 = t155 * t171 + t156 * t158;
t149 = pkin(11) + qJ(4);
t146 = sin(t149);
t147 = cos(t149);
t151 = sin(pkin(6));
t176 = t151 * t196;
t118 = t133 * t147 - t146 * t176;
t132 = t155 * t156 - t158 * t171;
t154 = sin(qJ(5));
t157 = cos(qJ(5));
t104 = t118 * t154 - t132 * t157;
t105 = t118 * t157 + t132 * t154;
t174 = t156 * t189;
t134 = t196 * t155 + t158 * t174;
t199 = -g(1) * t134 - g(2) * t132;
t190 = g(3) * t151;
t186 = t147 * t154;
t185 = t147 * t157;
t184 = t151 * t155;
t183 = t151 * t156;
t182 = t151 * t158;
t181 = t157 * t158;
t180 = t196 * pkin(1) + pkin(8) * t183;
t150 = sin(pkin(11));
t178 = t150 * t183;
t177 = t154 * t182;
t175 = -t156 * pkin(1) + pkin(8) * t176;
t173 = -t133 * t146 - t147 * t176;
t172 = t150 * t176;
t135 = -t155 * t174 + t196 * t158;
t122 = t135 * t147 + t146 * t183;
t108 = t122 * t154 - t134 * t157;
t126 = t189 * t146 + t147 * t184;
t115 = t126 * t154 + t151 * t181;
t165 = g(1) * t108 + g(2) * t104 + g(3) * t115;
t160 = g(3) * t182 + t199;
t153 = -pkin(9) - qJ(3);
t152 = cos(pkin(11));
t145 = pkin(3) * t152 + pkin(2);
t124 = (t147 * t181 + t154 * t155) * t151;
t123 = t147 * t177 - t157 * t184;
t121 = t135 * t146 - t147 * t183;
t116 = t126 * t157 - t177;
t113 = -t134 * t185 + t135 * t154;
t112 = -t134 * t186 - t135 * t157;
t111 = -t132 * t185 + t133 * t154;
t110 = -t132 * t186 - t133 * t157;
t109 = t122 * t157 + t134 * t154;
t1 = [(g(1) * t156 - g(2) * t196) * MDP(2) + (g(1) * t196 + g(2) * t156) * MDP(3) + (g(1) * t133 - g(2) * t135) * MDP(9) + (-g(1) * (-t133 * t152 + t172) - g(2) * (t135 * t152 + t178)) * MDP(11) + (-g(1) * (t133 * t150 + t152 * t176) - g(2) * (-t135 * t150 + t152 * t183)) * MDP(12) + (-g(1) * (-pkin(2) * t133 - qJ(3) * t132 + t175) - g(2) * (pkin(2) * t135 + qJ(3) * t134 + t180)) * MDP(14) + (g(1) * t118 - g(2) * t122) * MDP(20) + (-g(1) * (pkin(3) * t172 - pkin(4) * t118 - pkin(5) * t105 + pkin(10) * t173 - qJ(6) * t104 + t132 * t153 - t133 * t145 + t175) - g(2) * (pkin(3) * t178 + pkin(4) * t122 + pkin(5) * t109 + pkin(10) * t121 + qJ(6) * t108 - t134 * t153 + t135 * t145 + t180)) * MDP(32) + t179 * (g(1) * t105 - g(2) * t109) + t197 * (-g(1) * t104 + g(2) * t108) + t198 * (g(1) * t173 + g(2) * t121) - t200 * (g(1) * t132 - g(2) * t134); (-g(1) * (-pkin(2) * t134 + qJ(3) * t135) - g(2) * (-pkin(2) * t132 + qJ(3) * t133) - (pkin(2) * t158 + qJ(3) * t155) * t190) * MDP(14) + (-g(1) * (pkin(5) * t113 + qJ(6) * t112 - t135 * t153) - g(2) * (pkin(5) * t111 + qJ(6) * t110 - t133 * t153) - g(3) * (t124 * pkin(5) + t123 * qJ(6)) + t153 * t155 * t190 + (-t158 * t190 - t199) * (pkin(4) * t147 + pkin(10) * t146 + t145)) * MDP(32) + t197 * (g(1) * t112 + g(2) * t110 + g(3) * t123) + t200 * (g(1) * t135 + g(2) * t133 + g(3) * t184) + t179 * (-g(1) * t113 - g(2) * t111 - g(3) * t124) + (-t152 * MDP(11) + MDP(12) * t150 - t147 * MDP(20) + t146 * t198 - MDP(9)) * t160; (MDP(14) + MDP(32)) * t160; (-pkin(10) * MDP(32) + t198) * (g(1) * t122 + g(2) * t118 + g(3) * t126) + (MDP(20) + MDP(32) * (pkin(5) * t157 + qJ(6) * t154 + pkin(4)) + t179 * t157 - t197 * t154) * (-g(3) * (-t146 * t184 + t189 * t147) - g(2) * t173 + g(1) * t121); (-g(1) * (-pkin(5) * t108 + qJ(6) * t109) - g(2) * (-pkin(5) * t104 + qJ(6) * t105) - g(3) * (-pkin(5) * t115 + qJ(6) * t116)) * MDP(32) + t179 * t165 + t197 * (g(1) * t109 + g(2) * t105 + g(3) * t116); -t165 * MDP(32);];
taug  = t1;
