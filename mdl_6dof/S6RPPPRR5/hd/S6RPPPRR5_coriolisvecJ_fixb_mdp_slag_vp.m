% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPPRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:45
% EndTime: 2019-03-09 01:37:49
% DurationCPUTime: 1.55s
% Computational Cost: add. (915->233), mult. (1691->350), div. (0->0), fcn. (944->6), ass. (0->109)
t173 = sin(pkin(9));
t174 = cos(pkin(9));
t158 = qJD(2) * t174 + qJD(3) * t173;
t149 = qJD(1) * t158;
t256 = 0.2e1 * t149;
t178 = sin(qJ(5));
t180 = cos(qJ(5));
t255 = t178 * MDP(18) + t180 * MDP(19);
t179 = cos(qJ(6));
t177 = sin(qJ(6));
t223 = t177 * qJD(5);
t232 = qJD(1) * t178;
t153 = t179 * t232 + t223;
t254 = t153 * MDP(20);
t159 = qJD(2) * t173 - qJD(3) * t174;
t150 = t159 * qJD(1);
t171 = t178 ^ 2;
t253 = (-t180 ^ 2 + t171) * MDP(14);
t252 = MDP(6) * qJ(2) + MDP(5) + MDP(7);
t176 = -pkin(1) - qJ(3);
t168 = qJ(2) * qJD(1) + qJD(3);
t161 = pkin(3) * qJD(1) + t168;
t162 = t176 * qJD(1) + qJD(2);
t141 = t173 * t161 + t174 * t162;
t138 = qJD(1) * pkin(7) + t141;
t135 = qJD(4) * t178 + t138 * t180;
t130 = qJD(5) * t135 + t150 * t178;
t249 = t130 * t177;
t248 = t130 * t179;
t227 = qJD(6) * t177;
t210 = t178 * t227;
t217 = qJD(1) * qJD(5);
t208 = t180 * t217;
t216 = qJD(5) * qJD(6);
t237 = (t208 + t216) * t179;
t142 = -qJD(1) * t210 + t237;
t247 = t142 * t177;
t226 = qJD(6) * t179;
t209 = t178 * t226;
t143 = t177 * t216 + (t180 * t223 + t209) * qJD(1);
t246 = t143 * t180;
t213 = t177 * t232;
t221 = t179 * qJD(5);
t151 = t213 - t221;
t231 = qJD(1) * t180;
t163 = -qJD(6) + t231;
t245 = t151 * t163;
t244 = t151 * t178;
t243 = t153 * t163;
t242 = t177 * t163;
t241 = t177 * t180;
t240 = t178 * t142;
t239 = t179 * t163;
t238 = t179 * t180;
t175 = pkin(3) + qJ(2);
t236 = t173 * t175 + t174 * t176;
t181 = qJD(5) ^ 2;
t182 = qJD(1) ^ 2;
t234 = t181 + t182;
t233 = MDP(9) * qJD(1);
t148 = pkin(7) + t236;
t230 = qJD(5) * t148;
t229 = qJD(5) * t178;
t228 = qJD(5) * t180;
t225 = qJD(6) * t180;
t134 = qJD(4) * t180 - t138 * t178;
t131 = -qJD(5) * pkin(5) - t134;
t224 = t131 * qJD(6);
t219 = t182 * MDP(10);
t218 = t182 * MDP(11);
t214 = 0.2e1 * t217;
t212 = t163 * t227;
t211 = t163 * t226;
t207 = MDP(24) * t229;
t206 = t234 * t180;
t205 = t234 * t178;
t132 = qJD(5) * pkin(8) + t135;
t204 = t148 * t163 + t132;
t203 = -t142 * t180 + t153 * t229;
t140 = t161 * t174 - t173 * t162;
t202 = -t173 * t176 + t174 * t175;
t200 = -0.2e1 * t174 * t217;
t199 = 0.2e1 * t208;
t198 = t163 * t210;
t197 = t163 * t209;
t196 = pkin(5) * t178 - pkin(8) * t180;
t193 = -pkin(5) * t180 - pkin(8) * t178 - pkin(4);
t133 = t193 * qJD(1) - t140;
t128 = t132 * t179 + t133 * t177;
t195 = t132 * t177 - t133 * t179;
t194 = qJD(1) * t171 - t163 * t180;
t192 = t173 * t177 + t174 * t238;
t191 = -t173 * t179 + t174 * t241;
t190 = t173 * t238 - t174 * t177;
t189 = -t173 * t241 - t174 * t179;
t188 = -t148 * t181 + t256;
t137 = -qJD(1) * pkin(4) - t140;
t187 = qJD(5) * (qJD(1) * (-pkin(4) - t202) + t137 - t159);
t186 = t194 * t177;
t129 = qJD(5) * t134 + t150 * t180;
t185 = t131 * qJD(5) + qJD(6) * t133 + t159 * t163 + t129;
t145 = t196 * qJD(5) - t158;
t184 = (-t177 * t225 - t178 * t221) * t163 + t153 * t228 + t240;
t183 = -(t178 * t223 - t179 * t225) * t163 + t151 * t228 + t178 * t143;
t157 = t196 * qJD(1);
t144 = t193 - t202;
t139 = t145 * qJD(1);
t136 = t179 * t139;
t1 = [(qJD(2) * t168 - qJD(3) * t162 + (qJ(2) * qJD(2) - qJD(3) * t176) * qJD(1)) * MDP(9) + MDP(10) * t256 - 0.2e1 * MDP(11) * t150 + (t140 * t158 + t141 * t159 + t149 * t202 + t150 * t236) * MDP(12) - t214 * t253 + (-t210 * t153 + t179 * t240) * MDP(20) + (-t151 * t179 - t153 * t177) * t228 * MDP(21) + (t194 * t221 + t198 + t203) * MDP(22) + (t197 + t246 + (-t186 - t244) * qJD(5)) * MDP(23) - t231 * t207 + (-t207 - (-t144 * t227 + t145 * t179) * MDP(25) + (t144 * t226 + t145 * t177) * MDP(26)) * t163 + (t181 * MDP(15) + t188 * MDP(18) + t187 * MDP(19) + t221 * t254 + (t151 * t230 + t185 * t177 + t204 * t226 - t136) * MDP(25) + (t153 * t230 + (-t204 * qJD(6) + t139) * t177 + t185 * t179) * MDP(26)) * t180 + (MDP(13) * t199 - t181 * MDP(16) + t187 * MDP(18) - t188 * MDP(19) + (-t247 - t143 * t179 + (t151 * t177 - t153 * t179) * qJD(6)) * MDP(21) + (t179 * t224 + t249 + t148 * t143 + t159 * t151 + (-t148 * t242 + (t144 * t179 - t148 * t241) * qJD(1) - t195) * qJD(5)) * MDP(25) + (-t177 * t224 + t248 + t148 * t142 + t159 * t153 + (-t148 * t239 - (t144 * t177 + t148 * t238) * qJD(1) - t128) * qJD(5)) * MDP(26)) * t178 + 0.2e1 * (qJD(3) * MDP(8) + t252 * qJD(2)) * qJD(1); (-qJD(3) - t168) * t233 - t174 * t219 + t173 * t218 + (-t149 * t173 + t150 * t174 + (-t140 * t174 - t141 * t173) * qJD(1)) * MDP(12) + (t173 * t178 * t214 - t174 * t206) * MDP(18) + (t173 * t199 + t174 * t205) * MDP(19) + (t173 * t212 + t183 * t174 + (t189 * t163 + (-qJD(5) * t191 - t173 * t151) * t178) * qJD(1)) * MDP(25) + (t173 * t211 + t184 * t174 + (-t190 * t163 + (-qJD(5) * t192 - t173 * t153) * t178) * qJD(1)) * MDP(26) - t252 * t182; -t182 * MDP(8) + (qJD(2) + t162) * t233 - t173 * t219 - t174 * t218 + (t149 * t174 + t150 * t173 + (-t140 * t173 + t141 * t174) * qJD(1)) * MDP(12) + (-t173 * t206 + t178 * t200) * MDP(18) + (t173 * t205 + t180 * t200) * MDP(19) + (-t174 * t212 + t183 * t173 + (t191 * t163 + (qJD(5) * t189 + t174 * t151) * t178) * qJD(1)) * MDP(25) + (-t174 * t211 + t184 * t173 + (t192 * t163 + (-qJD(5) * t190 + t174 * t153) * t178) * qJD(1)) * MDP(26); (t197 - t246) * MDP(25) + (-t198 + t203) * MDP(26) - t255 * t181 + ((-t186 + t244) * MDP(25) - t194 * MDP(26) * t179) * qJD(5); (-t153 * t239 + t247) * MDP(20) + ((t142 + t245) * t179 + (-t143 + t243) * t177) * MDP(21) + (-t211 + (t163 * t238 + (-t153 + t223) * t178) * qJD(1)) * MDP(22) + (t212 + (-t163 * t241 + (t151 + t221) * t178) * qJD(1)) * MDP(23) + t163 * MDP(24) * t232 + (-pkin(5) * t143 - t248 + (-t134 * t177 + t157 * t179) * t163 - t135 * t151 + (pkin(8) * t239 + t131 * t177) * qJD(6) + (t195 * t178 + (-pkin(8) * t229 - t131 * t180) * t177) * qJD(1)) * MDP(25) + (-pkin(5) * t142 + t249 - (t134 * t179 + t157 * t177) * t163 - t135 * t153 + (-pkin(8) * t242 + t131 * t179) * qJD(6) + (-t131 * t238 + (-pkin(8) * t221 + t128) * t178) * qJD(1)) * MDP(26) + t255 * (-qJD(1) * t137 - t150) + (-t178 * t180 * MDP(13) + t253) * t182; t151 * t254 + (-t151 ^ 2 + t153 ^ 2) * MDP(21) + (t237 - t245) * MDP(22) + (-t177 * t208 - t243) * MDP(23) + qJD(1) * t207 + (-t128 * t163 - t177 * t129 - t131 * t153 + t136) * MDP(25) + (-t179 * t129 + t131 * t151 - t177 * t139 + t163 * t195) * MDP(26) + (-MDP(22) * t213 - t153 * MDP(23) - t128 * MDP(25) + t195 * MDP(26)) * qJD(6);];
tauc  = t1;
