% Calculate Gravitation load on the joints for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:39:42
% EndTime: 2019-03-09 08:39:44
% DurationCPUTime: 0.78s
% Computational Cost: add. (306->114), mult. (749->161), div. (0->0), fcn. (811->8), ass. (0->58)
t195 = MDP(11) + MDP(15);
t194 = -MDP(12) + MDP(17);
t190 = MDP(24) + MDP(26);
t189 = MDP(25) - MDP(28);
t193 = MDP(10) - MDP(13) - MDP(16) + MDP(27);
t149 = sin(qJ(2));
t150 = sin(qJ(1));
t153 = cos(qJ(1));
t166 = g(1) * t153 + g(2) * t150;
t192 = t149 * t166;
t146 = sin(pkin(9));
t147 = cos(pkin(9));
t148 = sin(qJ(5));
t151 = cos(qJ(5));
t163 = t146 * t148 + t147 * t151;
t191 = t163 * t149;
t188 = t149 * (-t146 * t151 + t147 * t148);
t187 = g(1) * t150;
t184 = g(3) * t149;
t152 = cos(qJ(2));
t142 = t152 * pkin(2);
t182 = t146 * t152;
t180 = t147 * t152;
t140 = t149 * qJ(3);
t179 = t149 * t150;
t178 = t149 * t153;
t177 = t150 * t152;
t176 = t152 * t153;
t175 = t153 * t146;
t174 = t142 + t140;
t173 = MDP(18) + MDP(29);
t172 = -pkin(1) - t142;
t171 = -qJ(4) * t146 - pkin(2);
t170 = pkin(3) * t180 + qJ(4) * t182 + t174;
t169 = t153 * pkin(1) + pkin(2) * t176 + t150 * pkin(7) + qJ(3) * t178;
t124 = t146 * t177 + t147 * t153;
t125 = t147 * t177 - t175;
t143 = t153 * pkin(7);
t164 = -t125 * pkin(3) - qJ(4) * t124 + t143;
t102 = t124 * t151 - t125 * t148;
t101 = t124 * t148 + t125 * t151;
t126 = -t150 * t147 + t152 * t175;
t127 = t150 * t146 + t147 * t176;
t104 = -t126 * t151 + t127 * t148;
t158 = g(1) * t104 - g(2) * t102 + g(3) * t188;
t155 = t127 * pkin(3) + t126 * qJ(4) + t169;
t154 = (t172 - t140) * t187;
t116 = -g(3) * t152 + t192;
t135 = qJ(3) * t176;
t132 = qJ(3) * t177;
t121 = t163 * t152;
t120 = t148 * t180 - t151 * t182;
t112 = t153 * t191;
t111 = t153 * t188;
t110 = t150 * t191;
t109 = t150 * t188;
t105 = t126 * t148 + t127 * t151;
t1 = [t166 * MDP(3) + (-g(1) * t143 - g(2) * t169 - t154) * MDP(14) + (-g(1) * t164 - g(2) * t155 - t154) * MDP(18) + (-g(1) * (-pkin(4) * t125 - pkin(5) * t101 + qJ(6) * t102 + t164) - g(2) * (t127 * pkin(4) + t105 * pkin(5) - pkin(8) * t178 + t104 * qJ(6) + t155) - ((pkin(8) - qJ(3)) * t149 + t172) * t187) * MDP(29) + t190 * (g(1) * t101 - g(2) * t105) + t189 * (g(1) * t102 + g(2) * t104) + t194 * (g(1) * t124 - g(2) * t126) + (t152 * MDP(9) + MDP(2)) * (-g(2) * t153 + t187) + t195 * (g(1) * t125 - g(2) * t127) - t193 * (g(1) * t179 - g(2) * t178); (-g(1) * (-pkin(2) * t178 + t135) - g(2) * (-pkin(2) * t179 + t132) - g(3) * t174) * MDP(14) + (-g(1) * t135 - g(2) * t132 - g(3) * t170 + (pkin(3) * t147 - t171) * t192) * MDP(18) + (-g(1) * (-t112 * pkin(5) - pkin(8) * t176 - t111 * qJ(6) + t135) - g(2) * (-pkin(5) * t110 - pkin(8) * t177 - qJ(6) * t109 + t132) - g(3) * (pkin(4) * t180 + pkin(5) * t121 + qJ(6) * t120 + t170) + (g(3) * pkin(8) + t166 * (-(-pkin(3) - pkin(4)) * t147 - t171)) * t149) * MDP(29) + t190 * (g(1) * t112 + g(2) * t110 - g(3) * t121) - t189 * (g(1) * t111 + g(2) * t109 - g(3) * t120) + t193 * (t152 * t166 + t184) + (t194 * t146 + t195 * t147 + MDP(9)) * t116; (-MDP(14) - t173) * t116; t173 * (-g(1) * t126 - g(2) * t124 - t146 * t184); (-g(1) * (-pkin(5) * t104 + qJ(6) * t105) - g(2) * (pkin(5) * t102 + qJ(6) * t101) - g(3) * (-pkin(5) * t188 + qJ(6) * t191)) * MDP(29) + t190 * t158 + t189 * (g(1) * t105 + g(2) * t101 + g(3) * t191); -t158 * MDP(29);];
taug  = t1;
