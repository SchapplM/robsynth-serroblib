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
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:57
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:54:51
% EndTime: 2021-01-16 02:54:55
% DurationCPUTime: 0.90s
% Computational Cost: add. (501->109), mult. (905->169), div. (0->0), fcn. (1067->12), ass. (0->62)
t156 = sin(qJ(3));
t159 = cos(qJ(3));
t153 = sin(pkin(6));
t157 = sin(qJ(2));
t187 = t153 * t157;
t194 = cos(pkin(6));
t207 = -t156 * t187 + t194 * t159;
t152 = sin(pkin(10));
t178 = t152 * t194;
t193 = cos(pkin(10));
t200 = cos(qJ(2));
t134 = t157 * t178 - t193 * t200;
t186 = t153 * t159;
t206 = t134 * t156 + t152 * t186;
t205 = -MDP(13) + MDP(24);
t204 = MDP(21) + MDP(23);
t203 = MDP(22) - MDP(25);
t180 = t153 * t200;
t170 = t194 * t193;
t136 = t152 * t200 + t157 * t170;
t177 = t153 * t193;
t202 = -t136 * t156 - t159 * t177;
t151 = qJ(3) + pkin(11);
t149 = sin(t151);
t150 = cos(t151);
t201 = pkin(4) * t150 + pkin(9) * t149;
t155 = sin(qJ(5));
t190 = t150 * t155;
t158 = cos(qJ(5));
t189 = t150 * t158;
t188 = t152 * t153;
t137 = t193 * t157 + t178 * t200;
t148 = pkin(3) * t159 + pkin(2);
t154 = qJ(4) + pkin(8);
t185 = -t134 * t154 - t137 * t148;
t184 = t148 * t180 + t154 * t187;
t183 = MDP(15) + MDP(26);
t179 = t158 * t200;
t135 = t152 * t157 - t170 * t200;
t175 = -t135 * t148 + t136 * t154;
t116 = t134 * t150 - t149 * t188;
t174 = t155 * t180;
t173 = t206 * pkin(3);
t169 = t207 * pkin(3);
t118 = t136 * t150 - t149 * t177;
t107 = t118 * t155 - t135 * t158;
t109 = -t116 * t155 - t137 * t158;
t128 = t149 * t194 + t150 * t187;
t121 = t128 * t155 + t153 * t179;
t167 = g(1) * t109 + g(2) * t107 + g(3) * t121;
t161 = g(2) * t202;
t160 = g(1) * t137 + g(2) * t135 - g(3) * t180;
t124 = (t150 * t179 + t155 * t157) * t153;
t123 = t150 * t174 - t158 * t187;
t122 = t128 * t158 - t174;
t114 = -t134 * t155 - t137 * t189;
t113 = t134 * t158 - t137 * t190;
t112 = -t135 * t189 + t136 * t155;
t111 = -t135 * t190 - t136 * t158;
t110 = -t116 * t158 + t137 * t155;
t108 = t118 * t158 + t135 * t155;
t1 = [(-MDP(1) - t183) * g(3); (-g(1) * t185 - g(2) * t175 - g(3) * t184) * MDP(15) + (-g(1) * (pkin(5) * t114 + qJ(6) * t113 - t137 * t201 + t185) - g(2) * (pkin(5) * t112 + qJ(6) * t111 - t135 * t201 + t175) - g(3) * (t124 * pkin(5) + t123 * qJ(6) + t201 * t180 + t184)) * MDP(26) + t203 * (g(1) * t113 + g(2) * t111 + g(3) * t123) + (MDP(4) - MDP(14)) * (-g(1) * t134 + g(2) * t136 + g(3) * t187) + t204 * (-g(1) * t114 - g(2) * t112 - g(3) * t124) + (MDP(10) * t159 - MDP(11) * t156 + MDP(12) * t150 + t149 * t205 + MDP(3)) * t160; (-g(1) * t206 - g(3) * t207 - t161) * MDP(10) + (-g(1) * (t134 * t159 - t156 * t188) - g(2) * (-t136 * t159 + t156 * t177) - g(3) * (-t156 * t194 - t157 * t186)) * MDP(11) + (-pkin(3) * t161 - g(1) * t173 - g(3) * t169) * MDP(15) + (-g(1) * (-pkin(9) * t116 + t173) - g(2) * (t202 * pkin(3) + t118 * pkin(9)) - g(3) * (pkin(9) * t128 + t169)) * MDP(26) + t205 * (g(1) * t116 - g(2) * t118 - g(3) * t128) + (-MDP(12) + (-pkin(5) * t158 - qJ(6) * t155 - pkin(4)) * MDP(26) - t204 * t158 + t203 * t155) * (g(3) * (-t149 * t187 + t150 * t194) + g(2) * (-t136 * t149 - t150 * t177) + g(1) * (t134 * t149 + t150 * t188)); -t183 * t160; (-g(1) * (-pkin(5) * t109 + qJ(6) * t110) - g(2) * (-pkin(5) * t107 + qJ(6) * t108) - g(3) * (-pkin(5) * t121 + qJ(6) * t122)) * MDP(26) + t204 * t167 + t203 * (g(1) * t110 + g(2) * t108 + g(3) * t122); -t167 * MDP(26);];
taug = t1;
