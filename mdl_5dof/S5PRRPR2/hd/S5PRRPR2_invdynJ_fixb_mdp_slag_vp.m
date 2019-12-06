% Calculate vector of inverse dynamics joint torques for
% S5PRRPR2
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:41
% EndTime: 2019-12-05 16:17:44
% DurationCPUTime: 1.20s
% Computational Cost: add. (892->211), mult. (1221->290), div. (0->0), fcn. (709->10), ass. (0->107)
t196 = sin(qJ(3));
t234 = qJD(3) * t196;
t224 = pkin(2) * t234;
t198 = cos(qJ(3));
t262 = pkin(2) * t198;
t240 = -qJD(2) * t224 + qJDD(2) * t262;
t219 = qJDD(4) - t240;
t186 = qJDD(2) + qJDD(3);
t261 = pkin(3) * t186;
t156 = t219 - t261;
t189 = pkin(8) + qJ(2);
t184 = qJ(3) + t189;
t175 = cos(t184);
t260 = g(2) * t175;
t264 = t156 + t260;
t190 = qJD(2) + qJD(3);
t228 = qJDD(2) * t196;
t233 = qJD(3) * t198;
t254 = t186 * qJ(4);
t148 = t254 + t190 * qJD(4) + (qJD(2) * t233 + t228) * pkin(2);
t193 = sin(pkin(9));
t194 = cos(pkin(9));
t144 = qJDD(1) * t193 + t148 * t194;
t143 = -t194 * qJDD(1) + t148 * t193;
t258 = t143 * t193;
t210 = t144 * t194 + t258;
t174 = sin(t184);
t171 = g(1) * t174;
t263 = -t171 - t240;
t195 = sin(qJ(5));
t197 = cos(qJ(5));
t205 = MDP(17) * t195 + MDP(18) * t197;
t235 = qJD(2) * t198;
t212 = -pkin(2) * t235 + qJD(4);
t259 = pkin(2) * qJD(2);
t226 = t196 * t259;
t160 = qJ(4) * t190 + t226;
t154 = -t194 * qJD(1) + t160 * t193;
t257 = t154 * t193;
t168 = t233 * pkin(2) + qJD(4);
t256 = t168 * t190;
t255 = t168 * t195;
t253 = t186 * t194;
t252 = t186 * t195;
t251 = t190 * t194;
t250 = t193 * t195;
t249 = t194 * t195;
t248 = t194 * t197;
t247 = t194 * t198;
t164 = -qJDD(5) + t253;
t246 = t195 * t164;
t245 = t195 * t197;
t244 = t197 * t164;
t243 = t264 * t193;
t242 = t175 * pkin(3) + t174 * qJ(4);
t241 = g(1) * t175 + g(2) * t174;
t187 = t193 ^ 2;
t188 = t194 ^ 2;
t239 = t187 + t188;
t192 = t197 ^ 2;
t238 = t195 ^ 2 - t192;
t232 = qJD(5) * t195;
t231 = qJD(5) * t197;
t230 = t164 * MDP(16);
t165 = -qJD(5) + t251;
t229 = -qJD(5) - t165;
t227 = t165 * t259;
t223 = qJ(4) * t231;
t222 = t190 * t234;
t221 = t190 * t231;
t220 = t165 * t232;
t218 = -pkin(3) * t174 + t175 * qJ(4);
t217 = t239 * t186;
t161 = -t194 * pkin(4) - t193 * pkin(7) - pkin(3);
t142 = t161 * t186 + t219;
t216 = t197 * t142 - t195 * t144;
t215 = t164 - t253;
t214 = t164 + t253;
t213 = t190 * t226;
t211 = t195 * t142 + t197 * t144;
t145 = t161 * t190 + t212;
t155 = qJD(1) * t193 + t160 * t194;
t209 = t145 * t197 - t155 * t195;
t208 = -t145 * t195 - t155 * t197;
t207 = t155 * t194 + t257;
t206 = t260 + t263;
t204 = -t241 + t210;
t150 = t174 * t249 + t175 * t197;
t152 = t174 * t197 - t175 * t249;
t203 = -g(1) * t150 - g(2) * t152 + (qJD(5) * t209 + t211) * t194 + t197 * t258;
t151 = -t174 * t248 + t175 * t195;
t153 = t174 * t195 + t175 * t248;
t202 = -g(1) * t151 - g(2) * t153 + t143 * t250 + t231 * t257;
t201 = t213 - t260;
t177 = -pkin(3) - t262;
t200 = t222 * pkin(2) + t177 * t186;
t158 = t193 * t232 * t251;
t199 = t194 * t230 + (t158 + (-t197 * t214 + t220) * t193) * MDP(14) + (t214 * t195 + (t165 + t251) * t231) * t193 * MDP(15) + t186 * MDP(5) + (0.2e1 * (qJD(5) * t190 * t238 - t186 * t245) * MDP(13) + (t186 * t192 - 0.2e1 * t195 * t221) * MDP(12)) * t187;
t185 = t190 ^ 2;
t183 = cos(t189);
t182 = sin(t189);
t176 = pkin(2) * t196 + qJ(4);
t163 = t194 * t171;
t159 = -t190 * pkin(3) + t212;
t157 = t161 - t262;
t134 = qJD(5) * t208 + t216;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (-t143 * t194 - g(3)) * MDP(11) + t158 * MDP(18) + (t144 * MDP(11) + (-qJD(5) * t165 * MDP(18) + MDP(17) * t215) * t195 + (t215 * MDP(18) + (t165 - t251) * MDP(17) * qJD(5)) * t197) * t193; qJDD(2) * MDP(2) + (g(1) * t182 - g(2) * t183) * MDP(3) + (g(1) * t183 + g(2) * t182) * MDP(4) + ((t186 * t198 - t222) * pkin(2) - t206) * MDP(6) + (((-qJDD(2) - t186) * t196 + (-qJD(2) - t190) * t233) * pkin(2) + t241) * MDP(7) + (t163 + (-t200 - t264) * t194) * MDP(8) + ((t200 - t171) * t193 + t243) * MDP(9) + (t176 * t217 + t239 * t256 + t204) * MDP(10) + (t156 * t177 + t159 * t224 - g(1) * (-pkin(2) * t182 + t218) - g(2) * (pkin(2) * t183 + t242) + t210 * t176 + t207 * t168) * MDP(11) + (-(-t157 * t232 + t197 * t224) * t165 - t157 * t244 + (-(-t176 * t231 - t255) * t165 + t176 * t246 - t134) * t194 + (t190 * t255 + (t221 + t252) * t176) * t187 + t202) * MDP(17) + ((t168 * t248 + t195 * t224) * t165 + (t157 * t195 + t176 * t248) * t164 + (t176 * t186 + t256) * t197 * t187 + (t197 * t157 * t165 + (-t257 + (-t165 * t194 - t187 * t190) * t176) * t195) * qJD(5) + t203) * MDP(18) + t199; (t201 - t263) * MDP(6) + ((-t228 + (-qJD(3) + t190) * t235) * pkin(2) + t241) * MDP(7) + (t163 + (-t156 + t201 + t261) * t194) * MDP(8) + ((-t171 - t213 - t261) * t193 + t243) * MDP(9) + (t212 * t190 * t239 + qJ(4) * t217 + t204) * MDP(10) + (-t156 * pkin(3) - g(1) * t218 - g(2) * t242 + t207 * qJD(4) + t210 * qJ(4) + (-t159 * t196 - t198 * t207) * t259) * MDP(11) + ((t220 - t244) * t161 + (-(-qJD(4) * t195 - t223) * t165 + qJ(4) * t246 - t134) * t194 + (-t195 * t247 + t196 * t197) * t227 + (qJ(4) * t252 + (t195 * t212 + t223) * t190) * t187 + t202) * MDP(17) + (qJD(4) * t165 * t248 + (qJ(4) * t248 + t161 * t195) * t164 + ((-qJ(4) * t249 + t161 * t197) * t165 - t154 * t250) * qJD(5) - (t195 * t196 + t197 * t247) * t227 + (t197 * t254 + (-qJ(4) * t232 + t197 * t212) * t190) * t187 + t203) * MDP(18) + t199; (-t190 * t207 + qJDD(4) + t206) * MDP(11) + (-pkin(3) * MDP(11) - MDP(8) * t194 + MDP(9) * t193) * t186 + (-MDP(17) * t197 + MDP(18) * t195) * t164 + (-t188 * MDP(10) + (-MDP(10) - t205) * t187) * t185 - t205 * t165 ^ 2; -t230 + (-g(1) * t152 + g(2) * t150 + t165 * t208 + t216) * MDP(17) + (g(1) * t153 - g(2) * t151 - t165 * t209 - t211) * MDP(18) + (MDP(17) * t208 - MDP(18) * t209) * qJD(5) + (MDP(12) * t245 - MDP(13) * t238) * t187 * t185 + ((MDP(14) * t197 - MDP(15) * t195) * t186 + t205 * g(3) + ((MDP(15) * t229 - t154 * MDP(17)) * t197 + (MDP(14) * t229 + t154 * MDP(18)) * t195) * t190) * t193;];
tau = t1;
