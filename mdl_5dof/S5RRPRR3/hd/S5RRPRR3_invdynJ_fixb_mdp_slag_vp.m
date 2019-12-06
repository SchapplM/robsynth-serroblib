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
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:30:36
% EndTime: 2019-12-05 18:30:38
% DurationCPUTime: 0.85s
% Computational Cost: add. (1182->164), mult. (1950->215), div. (0->0), fcn. (1145->14), ass. (0->98)
t208 = qJ(1) + qJ(2);
t195 = pkin(9) + qJ(4) + t208;
t187 = sin(t195);
t188 = cos(t195);
t273 = g(2) * t188 + g(3) * t187;
t210 = cos(pkin(9));
t191 = pkin(2) * t210 + pkin(3);
t212 = sin(qJ(4));
t216 = cos(qJ(4));
t209 = sin(pkin(9));
t266 = pkin(2) * t209;
t232 = t191 * t216 - t212 * t266;
t167 = -pkin(4) - t232;
t253 = t212 * t191 + t216 * t266;
t168 = pkin(8) + t253;
t204 = qJDD(1) + qJDD(2);
t200 = qJDD(4) + t204;
t219 = qJD(5) ^ 2;
t217 = cos(qJ(2));
t213 = sin(qJ(2));
t259 = t210 * t213;
t230 = pkin(1) * (-t209 * t217 - t259);
t169 = qJD(1) * t230;
t260 = t209 * t213;
t229 = pkin(1) * (t210 * t217 - t260);
t171 = qJD(1) * t229;
t205 = qJD(1) + qJD(2);
t201 = qJD(4) + t205;
t243 = (-t253 * qJD(4) - t169 * t216 + t171 * t212) * t201;
t272 = -t167 * t200 - t168 * t219 + t243;
t248 = qJDD(1) * t213;
t251 = qJD(1) * t217;
t271 = pkin(1) * (qJD(2) * t251 + t248);
t267 = pkin(1) * t217;
t196 = pkin(2) + t267;
t237 = -pkin(1) * t260 + t210 * t196;
t166 = pkin(3) + t237;
t173 = pkin(1) * t259 + t196 * t209;
t256 = t212 * t166 + t216 * t173;
t202 = sin(t208);
t203 = cos(t208);
t270 = g(2) * t203 + g(3) * t202;
t179 = pkin(1) * t251 + pkin(2) * t205;
t268 = pkin(1) * t213;
t247 = qJD(1) * t268;
t157 = t210 * t179 - t209 * t247;
t155 = pkin(3) * t205 + t157;
t158 = t179 * t209 + t210 * t247;
t145 = t155 * t212 + t158 * t216;
t197 = qJDD(1) * t267;
t165 = pkin(2) * t204 - qJD(2) * t247 + t197;
t152 = t210 * t165 - t209 * t271;
t147 = pkin(3) * t204 + t152;
t153 = t165 * t209 + t210 * t271;
t269 = t145 * qJD(4) - t216 * t147 + t212 * t153;
t261 = t158 * t212;
t221 = g(3) * t188 - (qJD(4) * t155 + t153) * t216 + qJD(4) * t261 - t212 * t147 - g(2) * t187;
t265 = pkin(4) * t200;
t170 = qJD(2) * t230;
t172 = qJD(2) * t229;
t263 = (t256 * qJD(4) - t170 * t216 + t172 * t212) * t201;
t262 = t145 * t201;
t258 = qJDD(3) - g(1);
t137 = -t265 + t269;
t144 = t155 * t216 - t261;
t142 = -pkin(4) * t201 - t144;
t211 = sin(qJ(5));
t215 = cos(qJ(5));
t249 = qJD(5) * t215;
t257 = t137 * t211 + t142 * t249;
t255 = -t232 * qJD(4) + t169 * t212 + t171 * t216;
t206 = t211 ^ 2;
t252 = -t215 ^ 2 + t206;
t250 = qJD(5) * t201;
t246 = t142 * qJD(5) * t211 + t273 * t215;
t245 = t197 + t270;
t244 = -g(2) * t202 + g(3) * t203;
t239 = qJD(1) * (-qJD(2) + t205);
t238 = qJD(2) * (-qJD(1) - t205);
t180 = qJDD(5) * t211 + t215 * t219;
t181 = qJDD(5) * t215 - t211 * t219;
t235 = 0.2e1 * (t200 * t211 * t215 - t252 * t250) * MDP(12) + (0.2e1 * t201 * t211 * t249 + t200 * t206) * MDP(11) + t180 * MDP(13) + t181 * MDP(14) + t200 * MDP(8);
t234 = t166 * t216 - t173 * t212;
t231 = t204 * MDP(4) + t235;
t227 = -pkin(8) * t219 + t262 + t265;
t148 = -pkin(4) - t234;
t149 = pkin(8) + t256;
t226 = -t148 * t200 - t149 * t219 - t263;
t225 = -pkin(8) * t200 - t142 * t201 + t221;
t138 = t234 * qJD(4) + t170 * t212 + t172 * t216;
t224 = -qJDD(5) * t149 + (t148 * t201 - t138) * qJD(5);
t223 = -pkin(4) * t250 - pkin(8) * qJDD(5) + qJD(5) * t144;
t222 = -qJDD(5) * t168 + (t167 * t201 + t255) * qJD(5);
t220 = -t269 + t273;
t218 = cos(qJ(1));
t214 = sin(qJ(1));
t199 = t201 ^ 2;
t1 = [qJDD(1) * MDP(1) + (g(2) * t218 + g(3) * t214) * MDP(2) + (-g(2) * t214 + g(3) * t218) * MDP(3) + ((t204 * t217 + t213 * t238) * pkin(1) + t245) * MDP(5) + (((-qJDD(1) - t204) * t213 + t217 * t238) * pkin(1) + t244) * MDP(6) + (t153 * t173 + t158 * t172 + t152 * t237 + t157 * t170 - g(2) * (-pkin(1) * t218 - pkin(2) * t203) - g(3) * (-pkin(1) * t214 - pkin(2) * t202)) * MDP(7) + (t234 * t200 + t220 - t263) * MDP(9) + (-t138 * t201 - t256 * t200 + t221) * MDP(10) + (t224 * t211 + (-t137 + t226) * t215 + t246) * MDP(16) + (t224 * t215 + (-t226 - t273) * t211 + t257) * MDP(17) + t231; (t239 * t268 + t245) * MDP(5) + ((t217 * t239 - t248) * pkin(1) + t244) * MDP(6) + (-t157 * t169 - t158 * t171 + (t152 * t210 + t153 * t209 + t270) * pkin(2)) * MDP(7) + (t232 * t200 + t220 + t243) * MDP(9) + (-t253 * t200 + t255 * t201 + t221) * MDP(10) + (t222 * t211 + (-t137 + t272) * t215 + t246) * MDP(16) + (t222 * t215 + (-t273 - t272) * t211 + t257) * MDP(17) + t231; t181 * MDP(16) - t180 * MDP(17) + t258 * MDP(7); (t220 + t262) * MDP(9) + (t144 * t201 + t221) * MDP(10) + t246 * MDP(16) + t257 * MDP(17) + ((-t137 + t227) * MDP(16) + t223 * MDP(17)) * t215 + (t223 * MDP(16) + (-t227 - t273) * MDP(17)) * t211 + t235; qJDD(5) * MDP(15) + t252 * MDP(12) * t199 + (t200 * MDP(14) + t258 * MDP(16) + t225 * MDP(17)) * t215 + (-t199 * t215 * MDP(11) + t200 * MDP(13) + t225 * MDP(16) - t258 * MDP(17)) * t211;];
tau = t1;
