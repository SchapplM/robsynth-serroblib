% Calculate vector of inverse dynamics joint torques for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4PRRR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:53
% EndTime: 2019-12-31 16:36:55
% DurationCPUTime: 1.74s
% Computational Cost: add. (668->230), mult. (1629->346), div. (0->0), fcn. (1277->10), ass. (0->114)
t197 = cos(qJ(3));
t239 = qJD(2) * qJD(3);
t230 = t197 * t239;
t194 = sin(qJ(3));
t236 = qJDD(2) * t194;
t210 = -t230 - t236;
t277 = -qJD(3) * qJD(4) + t210;
t195 = sin(qJ(2));
t190 = sin(pkin(4));
t255 = qJD(1) * t190;
t173 = qJD(2) * pkin(6) + t195 * t255;
t192 = cos(pkin(4));
t254 = qJD(1) * t192;
t157 = t173 * t197 + t194 * t254;
t198 = cos(qJ(2));
t240 = qJD(1) * qJD(2);
t160 = qJDD(2) * pkin(6) + (qJDD(1) * t195 + t198 * t240) * t190;
t237 = qJDD(1) * t192;
t227 = t197 * t237;
t139 = -qJDD(3) * pkin(3) + qJD(3) * t157 + t160 * t194 - t227;
t251 = qJD(2) * t197;
t180 = -qJD(4) + t251;
t189 = sin(pkin(8));
t191 = cos(pkin(8));
t262 = t192 * t195;
t162 = t189 * t198 + t191 * t262;
t164 = -t189 * t262 + t191 * t198;
t265 = t190 * t195;
t165 = -t192 * t197 + t194 * t265;
t264 = t190 * t197;
t208 = g(1) * (-t164 * t194 + t189 * t264) + g(2) * (-t162 * t194 - t191 * t264) - g(3) * t165;
t221 = pkin(3) * t194 - pkin(7) * t197;
t276 = (pkin(7) * qJD(4) + t221 * qJD(2)) * t180 - t139 - t208;
t199 = qJD(3) ^ 2;
t231 = t195 * t240;
t263 = t190 * t198;
t217 = -qJDD(1) * t263 + t190 * t231;
t261 = t192 * t198;
t219 = g(1) * (t189 * t261 + t191 * t195) + g(2) * (t189 * t195 - t191 * t261);
t275 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t199 + (-g(3) * t198 + t231) * t190 - t217 + t219;
t155 = qJD(3) * pkin(7) + t157;
t274 = (pkin(6) * t180 + t155) * qJD(4) + t219;
t273 = qJD(2) * pkin(2);
t193 = sin(qJ(4));
t196 = cos(qJ(4));
t222 = t193 * qJDD(3) - t277 * t196;
t244 = qJD(4) * t194;
t229 = qJD(2) * t244;
t144 = -t193 * t229 + t222;
t271 = t144 * t193;
t241 = t196 * qJD(3);
t253 = qJD(2) * t194;
t168 = t193 * t253 - t241;
t270 = t168 * t180;
t249 = qJD(3) * t193;
t252 = qJD(2) * t196;
t170 = t194 * t252 + t249;
t269 = t170 * t180;
t268 = t170 * t196;
t267 = t180 * t197;
t266 = t190 * t194;
t260 = t193 * t195;
t259 = t194 * t198;
t258 = t197 * t198;
t257 = qJDD(1) - g(3);
t187 = t194 ^ 2;
t256 = -t197 ^ 2 + t187;
t250 = qJD(3) * t170;
t248 = qJD(3) * t194;
t247 = qJD(3) * t197;
t175 = -pkin(3) * t197 - pkin(7) * t194 - pkin(2);
t246 = qJD(4) * t175;
t245 = qJD(4) * t193;
t243 = qJD(4) * t196;
t242 = qJD(4) * t198;
t235 = t197 * qJDD(2);
t234 = t198 * t255;
t233 = qJD(2) * t263;
t232 = t180 * t241;
t228 = t194 * t237;
t174 = -t234 - t273;
t224 = -qJD(2) * t174 - t160;
t156 = -t173 * t194 + t197 * t254;
t158 = t175 * qJD(2) - t234;
t223 = qJDD(3) * pkin(7) + qJD(3) * t156 + qJD(4) * t158 + t160 * t197 + t228;
t220 = g(1) * (t164 * t197 + t189 * t266) + g(2) * (t162 * t197 - t191 * t266);
t218 = g(1) * t164 + g(2) * t162;
t209 = -t194 * t239 + t235;
t167 = qJDD(4) - t209;
t172 = t221 * qJD(3);
t216 = -t175 * t167 + t172 * t180;
t166 = t192 * t194 + t195 * t264;
t214 = t196 * t258 + t260;
t213 = -t193 * t258 + t195 * t196;
t212 = t167 * t193 - t180 * t243;
t211 = t167 * t196 + t180 * t245;
t207 = -qJD(4) * t155 - t219;
t205 = t220 - t223;
t154 = -qJD(3) * pkin(3) - t156;
t204 = -pkin(6) * t167 + qJD(3) * t154 + t223;
t202 = -pkin(7) * t167 + (-t154 - t156) * t180;
t201 = -pkin(6) * qJDD(3) + (t174 + t234 - t273) * qJD(3);
t200 = qJD(2) ^ 2;
t184 = t196 * qJDD(3);
t153 = t166 * t196 - t193 * t263;
t152 = -t166 * t193 - t196 * t263;
t151 = t166 * qJD(3) + t194 * t233;
t150 = -t165 * qJD(3) + t197 * t233;
t145 = t196 * t229 - t184 + (t236 + (qJD(4) + t251) * qJD(3)) * t193;
t143 = qJD(2) * t172 + t175 * qJDD(2) + t217;
t142 = t196 * t143;
t141 = t155 * t196 + t158 * t193;
t140 = -t155 * t193 + t158 * t196;
t1 = [t257 * MDP(1) + (-qJD(3) * t151 - qJDD(3) * t165) * MDP(10) + (-qJD(3) * t150 - qJDD(3) * t166) * MDP(11) + (-(-t150 * t193 - t166 * t243) * t180 + t152 * t167 + t151 * t168 + t165 * t145) * MDP(17) + ((t150 * t196 - t166 * t245) * t180 - t153 * t167 + t151 * t170 + t165 * t144) * MDP(18) + ((-qJDD(2) * MDP(4) + (-MDP(10) * t197 + MDP(11) * t194 - MDP(3)) * t200) * t195 + (t209 * MDP(10) + t210 * MDP(11) + qJDD(2) * MDP(3) - t200 * MDP(4)) * t198 + (-(t193 * t242 + t195 * t252) * MDP(17) + (qJD(2) * t260 - t196 * t242) * MDP(18)) * t180) * t190; qJDD(2) * MDP(2) + (t257 * t263 + t219) * MDP(3) + (-t257 * t265 + t218) * MDP(4) + (qJDD(2) * t187 + 0.2e1 * t194 * t230) * MDP(5) + 0.2e1 * (t194 * t235 - t256 * t239) * MDP(6) + (qJDD(3) * t194 + t197 * t199) * MDP(7) + (qJDD(3) * t197 - t194 * t199) * MDP(8) + (t201 * t194 + t275 * t197) * MDP(10) + (-t275 * t194 + t201 * t197) * MDP(11) + (t144 * t194 * t196 + (-t193 * t244 + t197 * t241) * t170) * MDP(12) + ((-t168 * t196 - t170 * t193) * t247 + (-t271 - t145 * t196 + (t168 * t193 - t268) * qJD(4)) * t194) * MDP(13) + ((-t144 - t232) * t197 + (t211 + t250) * t194) * MDP(14) + ((t180 * t249 + t145) * t197 + (-qJD(3) * t168 - t212) * t194) * MDP(15) + (-t167 * t197 - t180 * t248) * MDP(16) + (t140 * t248 - t142 * t197 + (t194 * t145 + t168 * t247) * pkin(6) + (t154 * t244 + t274 * t197 - t216) * t196 + (-(pkin(6) * t248 - t246) * t180 + t139 * t194 + t204 * t197 - t218) * t193 + (-g(3) * t214 + (-t168 * t259 + t213 * t180) * qJD(1)) * t190) * MDP(17) + (t216 * t193 + (t180 * t246 - t218) * t196 + (pkin(6) * t250 + t204 * t196 + (t143 - t274) * t193) * t197 + (-t154 * t245 - t141 * qJD(3) + t139 * t196 + (t144 - t232) * pkin(6)) * t194 + (-g(3) * t213 + (-t170 * t259 - t214 * t180) * qJD(1)) * t190) * MDP(18); MDP(7) * t236 + MDP(8) * t235 + qJDD(3) * MDP(9) + (t224 * t194 - t208 + t227) * MDP(10) + (g(3) * t166 + t224 * t197 + t220 - t228) * MDP(11) + (-t180 * t268 + t271) * MDP(12) + ((t144 + t270) * t196 + (-t145 + t269) * t193) * MDP(13) + ((-t170 * t194 + t196 * t267) * qJD(2) + t212) * MDP(14) + ((t168 * t194 - t193 * t267) * qJD(2) + t211) * MDP(15) + t180 * MDP(16) * t253 + (-pkin(3) * t145 - t140 * t253 - t157 * t168 + t202 * t193 + t276 * t196) * MDP(17) + (-pkin(3) * t144 + t141 * t253 - t157 * t170 - t276 * t193 + t202 * t196) * MDP(18) + (-t194 * t197 * MDP(5) + t256 * MDP(6)) * t200; t170 * t168 * MDP(12) + (-t168 ^ 2 + t170 ^ 2) * MDP(13) + (t222 - t270) * MDP(14) + (t184 - t269) * MDP(15) + t167 * MDP(16) + (-g(3) * t152 - t141 * t180 - t154 * t170 + t142) * MDP(17) + (g(3) * t153 - t140 * t180 + t154 * t168) * MDP(18) + (-MDP(15) * t229 + t207 * MDP(17) + t205 * MDP(18)) * t196 + (-MDP(14) * t229 + t277 * MDP(15) + t205 * MDP(17) + (-t143 - t207) * MDP(18)) * t193;];
tau = t1;
