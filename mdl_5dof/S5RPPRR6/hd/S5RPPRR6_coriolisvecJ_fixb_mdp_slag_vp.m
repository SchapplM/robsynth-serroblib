% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:09
% EndTime: 2019-12-31 17:58:13
% DurationCPUTime: 1.22s
% Computational Cost: add. (973->177), mult. (2416->248), div. (0->0), fcn. (1772->8), ass. (0->85)
t166 = sin(pkin(9));
t168 = cos(pkin(9));
t190 = qJD(1) * (t166 ^ 2 + t168 ^ 2);
t218 = MDP(7) * t190;
t173 = cos(qJ(4));
t199 = qJD(1) * t173;
t157 = t168 * t199;
t171 = sin(qJ(4));
t200 = qJD(1) * t171;
t192 = t166 * t200;
t145 = t157 - t192;
t143 = qJD(5) - t145;
t152 = t166 * t171 - t173 * t168;
t177 = t152 * qJD(3);
t217 = qJD(1) * t177;
t153 = t166 * t173 + t168 * t171;
t146 = t153 * qJD(1);
t158 = sin(pkin(8)) * pkin(1) + qJ(3);
t155 = t158 * qJD(1);
t162 = t168 * qJD(2);
t213 = pkin(6) * qJD(1);
t132 = t162 + (-t155 - t213) * t166;
t140 = t166 * qJD(2) + t168 * t155;
t133 = t168 * t213 + t140;
t117 = t132 * t171 + t133 * t173;
t178 = t153 * qJD(3);
t112 = qJD(1) * t178 + qJD(4) * t117;
t216 = (pkin(4) * t146 + t143 * pkin(7)) * t143 + t112;
t116 = t132 * t173 - t133 * t171;
t111 = qJD(4) * t116 - t217;
t114 = -qJD(4) * pkin(4) - t116;
t214 = pkin(6) + t158;
t149 = t214 * t166;
t150 = t214 * t168;
t183 = -t149 * t173 - t150 * t171;
t118 = t183 * qJD(4) - t177;
t154 = -cos(pkin(8)) * pkin(1) - pkin(3) * t168 - pkin(2);
t144 = t154 * qJD(1) + qJD(3);
t120 = -pkin(4) * t145 - pkin(7) * t146 + t144;
t126 = pkin(4) * t152 - pkin(7) * t153 + t154;
t128 = -t149 * t171 + t150 * t173;
t148 = t153 * qJD(4);
t142 = qJD(1) * t148;
t147 = t152 * qJD(4);
t215 = t112 * t153 - t114 * t147 - t128 * t142 - (qJD(5) * t126 + t118) * t143 - (qJD(5) * t120 + t111) * t152;
t170 = sin(qJ(5));
t197 = qJD(5) * t170;
t156 = qJD(4) * t157;
t141 = -qJD(4) * t192 + t156;
t172 = cos(qJ(5));
t196 = t172 * qJD(4);
t202 = qJD(5) * t196 + t172 * t141;
t121 = -t146 * t197 + t202;
t212 = t121 * t170;
t211 = t126 * t142;
t204 = t146 * t170;
t134 = -t196 + t204;
t210 = t134 * t143;
t209 = t134 * t146;
t136 = qJD(4) * t170 + t146 * t172;
t208 = t136 * t143;
t207 = t136 * t146;
t206 = t141 * t170;
t205 = t142 * t170;
t138 = t172 * t142;
t203 = t121 * t152 + t136 * t148;
t198 = qJD(5) * t153;
t194 = t153 * t205;
t193 = t153 * t138;
t189 = t143 * t172;
t115 = qJD(4) * pkin(7) + t117;
t110 = t115 * t172 + t120 * t170;
t186 = t115 * t170 - t120 * t172;
t122 = qJD(5) * t136 + t206;
t185 = -t122 * t152 - t134 * t148;
t184 = (-t155 * t166 + t162) * t166 - t140 * t168;
t182 = t138 + (t145 * t170 - t197) * t143;
t181 = t147 * t170 - t172 * t198;
t180 = t147 * t172 + t153 * t197;
t175 = -pkin(7) * t142 + (t114 + t116) * t143;
t131 = pkin(4) * t148 + pkin(7) * t147;
t124 = pkin(4) * t142 - pkin(7) * t141;
t123 = t172 * t124;
t119 = t128 * qJD(4) + t178;
t1 = [(t141 * t153 - t146 * t147) * MDP(9) + (-t141 * t152 - t142 * t153 - t145 * t147 - t146 * t148) * MDP(10) + (t142 * t154 + t144 * t148) * MDP(14) + (t141 * t154 - t144 * t147) * MDP(15) + (t121 * t153 * t172 - t180 * t136) * MDP(16) + (-(-t134 * t172 - t136 * t170) * t147 + (-t212 - t122 * t172 + (t134 * t170 - t136 * t172) * qJD(5)) * t153) * MDP(17) + (-t180 * t143 + t193 + t203) * MDP(18) + (t181 * t143 + t185 - t194) * MDP(19) + (t142 * t152 + t143 * t148) * MDP(20) + (-t186 * t148 + t119 * t134 - t183 * t122 + t123 * t152 + (t211 + t131 * t143 + (t114 * t153 - t115 * t152 - t128 * t143) * qJD(5)) * t172 + t215 * t170) * MDP(21) + (-t110 * t148 + t119 * t136 - t183 * t121 + (-(-qJD(5) * t128 + t131) * t143 - t211 - (-qJD(5) * t115 + t124) * t152 - t114 * t198) * t170 + t215 * t172) * MDP(22) + (0.2e1 * t218 + (t158 * t190 - t184) * MDP(8)) * qJD(3) + (-t147 * MDP(11) - t148 * MDP(12) - t119 * MDP(14) - t118 * MDP(15)) * qJD(4); (-t185 - t194) * MDP(21) + (-t193 + t203) * MDP(22) + (t181 * MDP(21) + t180 * MDP(22)) * t143 + (-MDP(14) * t148 + MDP(15) * t147) * qJD(4); t156 * MDP(15) + (t182 - t209) * MDP(21) + (-t143 ^ 2 * t172 - t205 - t207) * MDP(22) + ((t166 * t199 + t168 * t200 + t146) * MDP(14) + (t145 - t192) * MDP(15)) * qJD(4) + (t184 * MDP(8) - t218) * qJD(1); -t145 ^ 2 * MDP(10) + (t156 + (-t145 - t192) * qJD(4)) * MDP(11) + (-t144 * t145 + t217) * MDP(15) + (t136 * t189 + t212) * MDP(16) + ((t121 - t210) * t172 + (-t122 - t208) * t170) * MDP(17) + (t143 * t189 + t205 - t207) * MDP(18) + (t182 + t209) * MDP(19) + (-pkin(4) * t122 - t117 * t134 + t175 * t170 - t172 * t216) * MDP(21) + (-pkin(4) * t121 - t117 * t136 + t170 * t216 + t175 * t172) * MDP(22) + (MDP(10) * t146 - t143 * MDP(20) + t186 * MDP(21) + t110 * MDP(22) - t145 * MDP(9) + (-qJD(3) - t144) * MDP(14)) * t146; t136 * t134 * MDP(16) + (-t134 ^ 2 + t136 ^ 2) * MDP(17) + (t202 + t210) * MDP(18) + (-t206 + t208) * MDP(19) + t142 * MDP(20) + (t110 * t143 - t111 * t170 - t114 * t136 + t123) * MDP(21) + (-t111 * t172 + t114 * t134 - t124 * t170 - t143 * t186) * MDP(22) + (-MDP(18) * t204 - t136 * MDP(19) - t110 * MDP(21) + t186 * MDP(22)) * qJD(5);];
tauc = t1;
