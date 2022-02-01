% Calculate vector of inverse dynamics joint torques for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRPRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:34:23
% EndTime: 2022-01-20 10:34:26
% DurationCPUTime: 0.94s
% Computational Cost: add. (1182->164), mult. (1950->215), div. (0->0), fcn. (1145->14), ass. (0->98)
t207 = qJD(1) + qJD(2);
t219 = cos(qJ(2));
t252 = qJD(1) * t219;
t181 = pkin(1) * t252 + pkin(2) * t207;
t211 = sin(pkin(9));
t212 = cos(pkin(9));
t215 = sin(qJ(2));
t271 = pkin(1) * t215;
t248 = qJD(1) * t271;
t159 = t212 * t181 - t211 * t248;
t157 = pkin(3) * t207 + t159;
t160 = t181 * t211 + t212 * t248;
t214 = sin(qJ(4));
t218 = cos(qJ(4));
t147 = t214 * t157 + t218 * t160;
t270 = pkin(1) * t219;
t199 = qJDD(1) * t270;
t206 = qJDD(1) + qJDD(2);
t167 = pkin(2) * t206 - qJD(2) * t248 + t199;
t249 = qJDD(1) * t215;
t274 = pkin(1) * (qJD(2) * t252 + t249);
t154 = t212 * t167 - t211 * t274;
t149 = pkin(3) * t206 + t154;
t155 = t167 * t211 + t212 * t274;
t210 = qJ(1) + qJ(2);
t197 = pkin(9) + qJ(4) + t210;
t190 = cos(t197);
t277 = -g(2) * t190 - t147 * qJD(4) + t218 * t149 - t214 * t155;
t202 = qJDD(4) + t206;
t268 = pkin(4) * t202;
t276 = -t268 - t277;
t193 = pkin(2) * t212 + pkin(3);
t269 = pkin(2) * t211;
t234 = t193 * t218 - t214 * t269;
t169 = -pkin(4) - t234;
t255 = t214 * t193 + t218 * t269;
t170 = pkin(8) + t255;
t221 = qJD(5) ^ 2;
t261 = t212 * t215;
t232 = pkin(1) * (-t211 * t219 - t261);
t171 = qJD(1) * t232;
t262 = t211 * t215;
t231 = pkin(1) * (t212 * t219 - t262);
t173 = qJD(1) * t231;
t203 = qJD(4) + t207;
t245 = (-t255 * qJD(4) - t171 * t218 + t173 * t214) * t203;
t275 = t169 * t202 + t170 * t221 - t245;
t198 = pkin(2) + t270;
t239 = -pkin(1) * t262 + t212 * t198;
t168 = pkin(3) + t239;
t175 = pkin(1) * t261 + t198 * t211;
t258 = t214 * t168 + t218 * t175;
t204 = sin(t210);
t205 = cos(t210);
t273 = g(1) * t204 - g(2) * t205;
t189 = sin(t197);
t263 = t160 * t214;
t223 = g(1) * t190 + g(2) * t189 - (qJD(4) * t157 + t155) * t218 + qJD(4) * t263 - t214 * t149;
t187 = g(1) * t189;
t172 = qJD(2) * t232;
t174 = qJD(2) * t231;
t265 = (t258 * qJD(4) - t172 * t218 + t174 * t214) * t203;
t264 = t147 * t203;
t260 = qJDD(3) - g(3);
t146 = t218 * t157 - t263;
t144 = -pkin(4) * t203 - t146;
t213 = sin(qJ(5));
t217 = cos(qJ(5));
t259 = t144 * qJD(5) * t213 + t217 * t187;
t257 = -t234 * qJD(4) + t171 * t214 + t173 * t218;
t254 = g(1) * t205 + g(2) * t204;
t208 = t213 ^ 2;
t253 = -t217 ^ 2 + t208;
t251 = qJD(5) * t203;
t250 = qJD(5) * t217;
t247 = t144 * t250 + t276 * t213;
t241 = qJD(1) * (-qJD(2) + t207);
t240 = qJD(2) * (-qJD(1) - t207);
t238 = t199 + t273;
t182 = qJDD(5) * t213 + t217 * t221;
t183 = qJDD(5) * t217 - t213 * t221;
t237 = 0.2e1 * (t202 * t213 * t217 - t253 * t251) * MDP(12) + (0.2e1 * t203 * t213 * t250 + t202 * t208) * MDP(11) + t182 * MDP(13) + t183 * MDP(14) + t202 * MDP(8);
t236 = t168 * t218 - t175 * t214;
t233 = t206 * MDP(4) + t237;
t229 = pkin(8) * t221 - t264 - t268;
t150 = -pkin(4) - t236;
t151 = pkin(8) + t258;
t228 = t150 * t202 + t151 * t221 + t265;
t227 = -pkin(8) * t202 - t144 * t203 + t223;
t140 = t236 * qJD(4) + t172 * t214 + t174 * t218;
t226 = -qJDD(5) * t151 + (t150 * t203 - t140) * qJD(5);
t225 = -pkin(4) * t251 - pkin(8) * qJDD(5) + qJD(5) * t146;
t224 = -qJDD(5) * t170 + (t169 * t203 + t257) * qJD(5);
t222 = t187 + t277;
t220 = cos(qJ(1));
t216 = sin(qJ(1));
t201 = t203 ^ 2;
t1 = [qJDD(1) * MDP(1) + (g(1) * t216 - g(2) * t220) * MDP(2) + (g(1) * t220 + g(2) * t216) * MDP(3) + ((t206 * t219 + t215 * t240) * pkin(1) + t238) * MDP(5) + (((-qJDD(1) - t206) * t215 + t219 * t240) * pkin(1) + t254) * MDP(6) + (t155 * t175 + t160 * t174 + t154 * t239 + t159 * t172 - g(1) * (-pkin(1) * t216 - pkin(2) * t204) - g(2) * (pkin(1) * t220 + pkin(2) * t205)) * MDP(7) + (t236 * t202 + t222 - t265) * MDP(9) + (-t140 * t203 - t258 * t202 + t223) * MDP(10) + (t226 * t213 + (-t228 - t276) * t217 + t259) * MDP(16) + (t226 * t217 + (t228 - t187) * t213 + t247) * MDP(17) + t233; (t241 * t271 + t238) * MDP(5) + ((t219 * t241 - t249) * pkin(1) + t254) * MDP(6) + (-t159 * t171 - t160 * t173 + (t154 * t212 + t155 * t211 + t273) * pkin(2)) * MDP(7) + (t234 * t202 + t222 + t245) * MDP(9) + (-t255 * t202 + t257 * t203 + t223) * MDP(10) + (t224 * t213 + (-t276 - t275) * t217 + t259) * MDP(16) + (t224 * t217 + (-t187 + t275) * t213 + t247) * MDP(17) + t233; t183 * MDP(16) - t182 * MDP(17) + t260 * MDP(7); (t222 + t264) * MDP(9) + (t146 * t203 + t223) * MDP(10) + t259 * MDP(16) + t247 * MDP(17) + (t225 * MDP(16) + (t229 - t187) * MDP(17)) * t213 + ((-t229 - t276) * MDP(16) + t225 * MDP(17)) * t217 + t237; qJDD(5) * MDP(15) + t253 * MDP(12) * t201 + (t202 * MDP(14) + t260 * MDP(16) + t227 * MDP(17)) * t217 + (-t201 * t217 * MDP(11) + t202 * MDP(13) + t227 * MDP(16) - t260 * MDP(17)) * t213;];
tau = t1;
