% Calculate Gravitation load on the joints for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:28:00
% EndTime: 2019-03-09 16:28:04
% DurationCPUTime: 1.14s
% Computational Cost: add. (468->142), mult. (1072->224), div. (0->0), fcn. (1267->12), ass. (0->59)
t194 = -MDP(16) + MDP(19) - MDP(24);
t192 = MDP(17) - MDP(20);
t149 = sin(qJ(3));
t197 = -qJ(4) * t149 - pkin(2);
t151 = sin(qJ(1));
t150 = sin(qJ(2));
t187 = cos(pkin(6));
t167 = t150 * t187;
t189 = cos(qJ(2));
t190 = cos(qJ(1));
t127 = t151 * t189 + t167 * t190;
t152 = cos(qJ(3));
t147 = sin(pkin(6));
t171 = t147 * t190;
t109 = t127 * t149 + t152 * t171;
t162 = t187 * t189;
t126 = t150 * t151 - t162 * t190;
t145 = pkin(11) + qJ(6);
t142 = sin(t145);
t143 = cos(t145);
t196 = t109 * t142 + t126 * t143;
t195 = t109 * t143 - t126 * t142;
t193 = MDP(10) - MDP(18);
t191 = pkin(4) + pkin(9);
t188 = g(3) * t147;
t183 = t126 * t152;
t128 = t150 * t190 + t151 * t162;
t182 = t128 * t152;
t181 = t142 * t149;
t180 = t143 * t149;
t146 = sin(pkin(11));
t179 = t146 * t149;
t178 = t147 * t150;
t177 = t147 * t151;
t176 = t147 * t152;
t148 = cos(pkin(11));
t175 = t148 * t149;
t173 = -pkin(3) * t183 + t197 * t126;
t172 = -pkin(3) * t182 + t197 * t128;
t170 = t147 * t189;
t169 = t149 * t189;
t168 = t152 * t189;
t110 = t127 * t152 - t149 * t171;
t166 = -t109 * pkin(3) + qJ(4) * t110;
t129 = -t151 * t167 + t189 * t190;
t113 = t129 * t149 - t151 * t176;
t114 = t129 * t152 + t149 * t177;
t165 = -t113 * pkin(3) + qJ(4) * t114;
t124 = t149 * t178 - t152 * t187;
t125 = t149 * t187 + t150 * t176;
t164 = -t124 * pkin(3) + qJ(4) * t125;
t163 = pkin(2) * t170 + pkin(9) * t178 + (pkin(3) * t168 + qJ(4) * t169) * t147;
t158 = t190 * pkin(1) + t129 * pkin(2) + t114 * pkin(3) + pkin(8) * t177 + qJ(4) * t113;
t157 = g(1) * t113 + g(2) * t109 + g(3) * t124;
t156 = g(1) * t114 + g(2) * t110 + g(3) * t125;
t154 = -pkin(1) * t151 - t127 * pkin(2) - pkin(3) * t110 + pkin(8) * t171 - qJ(4) * t109;
t102 = t113 * t142 + t128 * t143;
t101 = t113 * t143 - t128 * t142;
t1 = [(g(1) * t151 - g(2) * t190) * MDP(2) + (g(1) * t190 + g(2) * t151) * MDP(3) + (g(1) * t127 - g(2) * t129) * MDP(9) + (-g(1) * (-pkin(9) * t126 + t154) - g(2) * (pkin(9) * t128 + t158)) * MDP(21) + (-g(1) * (-t109 * t146 - t126 * t148) - g(2) * (t113 * t146 + t128 * t148)) * MDP(22) + (-g(1) * (-t109 * t148 + t126 * t146) - g(2) * (t113 * t148 - t128 * t146)) * MDP(23) + (-g(1) * (-qJ(5) * t110 - t126 * t191 + t154) - g(2) * (qJ(5) * t114 + t128 * t191 + t158)) * MDP(25) + (g(1) * t196 - g(2) * t102) * MDP(31) + (g(1) * t195 - g(2) * t101) * MDP(32) + t192 * (-g(1) * t109 + g(2) * t113) + t194 * (-g(1) * t110 + g(2) * t114) - t193 * (g(1) * t126 - g(2) * t128); (-g(1) * (t129 * pkin(9) + t172) - g(2) * (t127 * pkin(9) + t173) - g(3) * t163) * MDP(21) + (-g(1) * (-t128 * t179 + t129 * t148) - g(2) * (-t126 * t179 + t127 * t148) - (t146 * t169 + t148 * t150) * t188) * MDP(22) + (-g(1) * (-t128 * t175 - t129 * t146) - g(2) * (-t126 * t175 - t127 * t146) - (-t146 * t150 + t148 * t169) * t188) * MDP(23) + (-g(1) * (-qJ(5) * t182 + t129 * t191 + t172) - g(2) * (-qJ(5) * t183 + t127 * t191 + t173) - g(3) * ((pkin(4) * t150 + qJ(5) * t168) * t147 + t163)) * MDP(25) + (-g(1) * (-t128 * t181 + t129 * t143) - g(2) * (-t126 * t181 + t127 * t143) - (t142 * t169 + t143 * t150) * t188) * MDP(31) + (-g(1) * (-t128 * t180 - t129 * t142) - g(2) * (-t126 * t180 - t127 * t142) - (-t142 * t150 + t143 * t169) * t188) * MDP(32) + t193 * (g(1) * t129 + g(2) * t127 + g(3) * t178) + (t192 * t149 + t194 * t152 - MDP(9)) * (-g(1) * t128 - g(2) * t126 + g(3) * t170); (-g(1) * t165 - g(2) * t166 - g(3) * t164) * MDP(21) + (-g(1) * (-qJ(5) * t113 + t165) - g(2) * (-qJ(5) * t109 + t166) - g(3) * (-qJ(5) * t124 + t164)) * MDP(25) - t194 * t157 + (-MDP(22) * t146 - MDP(23) * t148 - MDP(31) * t142 - MDP(32) * t143 + t192) * t156; -(MDP(21) + MDP(25)) * t157; -t156 * MDP(25); (-g(1) * t101 - g(2) * t195 - g(3) * (t124 * t143 + t142 * t170)) * MDP(31) + (g(1) * t102 + g(2) * t196 - g(3) * (-t124 * t142 + t143 * t170)) * MDP(32);];
taug  = t1;
