% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRPRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:57:57
% EndTime: 2019-12-05 15:58:02
% DurationCPUTime: 1.68s
% Computational Cost: add. (1082->220), mult. (2934->322), div. (0->0), fcn. (2327->10), ass. (0->106)
t199 = cos(pkin(10));
t205 = cos(qJ(4));
t233 = qJD(2) * t205;
t190 = t199 * t233;
t197 = sin(pkin(10));
t202 = sin(qJ(4));
t235 = qJD(2) * t202;
t225 = t197 * t235;
t173 = -t190 + t225;
t170 = qJD(5) + t173;
t206 = cos(qJ(2));
t198 = sin(pkin(5));
t237 = qJD(1) * t198;
t226 = t206 * t237;
t179 = (qJD(3) + t226) * qJD(2);
t180 = t197 * t202 - t205 * t199;
t261 = t180 * t179;
t181 = t197 * t205 + t199 * t202;
t260 = t181 * t179;
t238 = t197 ^ 2 + t199 ^ 2;
t259 = t238 * MDP(7);
t218 = qJD(3) - t226;
t203 = sin(qJ(2));
t227 = t203 * t237;
t183 = qJD(2) * qJ(3) + t227;
t200 = cos(pkin(5));
t236 = qJD(1) * t200;
t189 = t199 * t236;
t255 = pkin(7) * qJD(2);
t154 = t189 + (-t183 - t255) * t197;
t164 = t199 * t183 + t197 * t236;
t155 = t199 * t255 + t164;
t140 = t154 * t202 + t155 * t205;
t134 = qJD(4) * t140 + t260;
t175 = qJD(2) * t181;
t258 = (pkin(4) * t175 + pkin(8) * t170) * t170 + t134;
t139 = t154 * t205 - t155 * t202;
t133 = qJD(4) * t139 - t261;
t135 = -qJD(4) * pkin(4) - t139;
t192 = -pkin(3) * t199 - pkin(2);
t169 = t192 * qJD(2) + t218;
t141 = pkin(4) * t173 - pkin(8) * t175 + t169;
t153 = pkin(4) * t180 - pkin(8) * t181 + t192;
t256 = pkin(7) + qJ(3);
t184 = t256 * t197;
t185 = t256 * t199;
t157 = -t184 * t202 + t185 * t205;
t177 = t181 * qJD(4);
t168 = qJD(2) * t177;
t176 = t180 * qJD(4);
t243 = t198 * t206;
t210 = t180 * t243;
t214 = -t184 * t205 - t185 * t202;
t241 = -qJD(1) * t210 + t180 * qJD(3) - t214 * qJD(4);
t257 = -(qJD(5) * t141 + t133) * t180 + t134 * t181 - t135 * t176 + (-qJD(5) * t153 + t241) * t170 - t157 * t168;
t254 = qJD(2) * pkin(2);
t253 = t135 * t181;
t201 = sin(qJ(5));
t232 = qJD(5) * t201;
t187 = qJD(4) * t190;
t167 = -qJD(4) * t225 + t187;
t204 = cos(qJ(5));
t229 = t204 * qJD(4);
t239 = qJD(5) * t229 + t204 * t167;
t144 = -t175 * t232 + t239;
t252 = t144 * t201;
t251 = t153 * t168;
t245 = t175 * t201;
t160 = -t229 + t245;
t250 = t160 * t170;
t249 = t160 * t175;
t162 = qJD(4) * t201 + t175 * t204;
t248 = t162 * t170;
t247 = t162 * t175;
t246 = t167 * t201;
t244 = t198 * t203;
t242 = t201 * t168;
t166 = t204 * t168;
t211 = t181 * t243;
t240 = -qJD(1) * t211 + t181 * qJD(3) + t157 * qJD(4);
t234 = qJD(2) * t203;
t231 = qJD(5) * t204;
t230 = qJD(5) * t206;
t223 = t238 * t179;
t222 = t170 * t204;
t219 = pkin(4) * t177 + pkin(8) * t176 - t227;
t136 = qJD(4) * pkin(8) + t140;
t132 = t136 * t204 + t141 * t201;
t217 = t136 * t201 - t141 * t204;
t216 = (-t183 * t197 + t189) * t197 - t164 * t199;
t171 = -t197 * t244 + t199 * t200;
t172 = t197 * t200 + t199 * t244;
t215 = t171 * t205 - t172 * t202;
t149 = t171 * t202 + t172 * t205;
t213 = t166 + (-t173 * t201 - t232) * t170;
t212 = -t176 * t204 - t181 * t232;
t209 = -pkin(8) * t168 + (t135 + t139) * t170;
t207 = qJD(2) ^ 2;
t186 = qJD(2) * t227;
t182 = t218 - t254;
t145 = t162 * qJD(5) + t246;
t143 = pkin(4) * t168 - pkin(8) * t167 + t186;
t142 = t204 * t143;
t138 = qJD(2) * t211 + t149 * qJD(4);
t137 = -qJD(2) * t210 + t215 * qJD(4);
t1 = [((-t137 * t201 - t149 * t231) * t170 - t149 * t242 + t138 * t160 - t215 * t145) * MDP(21) + (-(t137 * t204 - t149 * t232) * t170 - t149 * t166 + t138 * t162 - t215 * t144) * MDP(22) + (-t171 * t197 + t172 * t199) * MDP(8) * t179 + (-MDP(14) * t138 - MDP(15) * t137) * qJD(4) + ((-t168 * t206 + t173 * t234) * MDP(14) + (-t167 * t206 + t175 * t234) * MDP(15) + ((t201 * t230 + t204 * t234) * t170 - t206 * t166) * MDP(21) + (-(t201 * t234 - t204 * t230) * t170 + t206 * t242) * MDP(22) + (t182 * t203 + (-t216 - t227) * t206) * qJD(2) * MDP(8) + ((-MDP(4) + t259) * t206 + (-MDP(5) * t199 + MDP(6) * t197 - MDP(3)) * t203) * t207) * t198; (t218 * qJD(2) * t238 + t223) * MDP(7) + (-t216 * qJD(3) + qJ(3) * t223 + (t216 * t206 + (-t182 - t254) * t203) * t237) * MDP(8) + (t167 * t181 - t175 * t176) * MDP(9) + (-t167 * t180 - t168 * t181 + t173 * t176 - t175 * t177) * MDP(10) + (t168 * t192 + t169 * t177 + (qJD(2) * t180 - t173) * t227) * MDP(14) + (t167 * t192 - t169 * t176) * MDP(15) + (t144 * t181 * t204 + t212 * t162) * MDP(16) + (-(-t160 * t204 - t162 * t201) * t176 + (-t252 - t145 * t204 + (t160 * t201 - t162 * t204) * qJD(5)) * t181) * MDP(17) + (t144 * t180 + t162 * t177 + t181 * t166 + t212 * t170) * MDP(18) + (-t181 * t242 - t145 * t180 - t160 * t177 + (t176 * t201 - t181 * t231) * t170) * MDP(19) + (t168 * t180 + t170 * t177) * MDP(20) + (-t217 * t177 + t142 * t180 - t214 * t145 + t240 * t160 + (t251 + t219 * t170 + (-t136 * t180 - t157 * t170 + t253) * qJD(5)) * t204 + t257 * t201) * MDP(21) + (-t132 * t177 - t214 * t144 + t240 * t162 + (-t251 - (-qJD(5) * t136 + t143) * t180 - qJD(5) * t253 + (qJD(5) * t157 - t219) * t170) * t201 + t257 * t204) * MDP(22) + (-MDP(11) * t176 - t177 * MDP(12) - MDP(14) * t240 + MDP(15) * t241) * qJD(4); (t216 * qJD(2) + t186) * MDP(8) + t187 * MDP(15) + (t213 - t249) * MDP(21) + (-t170 ^ 2 * t204 - t242 - t247) * MDP(22) - t207 * t259 + ((t197 * t233 + t199 * t235 + t175) * MDP(14) + (-t173 - t225) * MDP(15)) * qJD(4); -t173 ^ 2 * MDP(10) + (t187 + (t173 - t225) * qJD(4)) * MDP(11) - MDP(14) * t260 + (t169 * t173 + t261) * MDP(15) + (t162 * t222 + t252) * MDP(16) + ((t144 - t250) * t204 + (-t145 - t248) * t201) * MDP(17) + (t170 * t222 + t242 - t247) * MDP(18) + (t213 + t249) * MDP(19) + (-pkin(4) * t145 - t140 * t160 + t209 * t201 - t258 * t204) * MDP(21) + (-pkin(4) * t144 - t140 * t162 + t258 * t201 + t209 * t204) * MDP(22) + (MDP(10) * t175 - t169 * MDP(14) - t170 * MDP(20) + MDP(21) * t217 + t132 * MDP(22) + t173 * MDP(9)) * t175; t162 * t160 * MDP(16) + (-t160 ^ 2 + t162 ^ 2) * MDP(17) + (t239 + t250) * MDP(18) + (-t246 + t248) * MDP(19) + t168 * MDP(20) + (t132 * t170 - t133 * t201 - t135 * t162 + t142) * MDP(21) + (-t133 * t204 + t135 * t160 - t143 * t201 - t170 * t217) * MDP(22) + (-MDP(18) * t245 - t162 * MDP(19) - t132 * MDP(21) + t217 * MDP(22)) * qJD(5);];
tauc = t1;
