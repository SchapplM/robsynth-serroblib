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
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:00:39
% EndTime: 2020-01-03 12:00:40
% DurationCPUTime: 0.83s
% Computational Cost: add. (1182->160), mult. (1950->213), div. (0->0), fcn. (1145->14), ass. (0->96)
t198 = qJD(1) + qJD(2);
t210 = cos(qJ(2));
t246 = qJD(1) * t210;
t176 = pkin(1) * t246 + pkin(2) * t198;
t202 = sin(pkin(9));
t203 = cos(pkin(9));
t206 = sin(qJ(2));
t264 = pkin(1) * t206;
t242 = qJD(1) * t264;
t154 = t203 * t176 - t202 * t242;
t152 = pkin(3) * t198 + t154;
t155 = t176 * t202 + t203 * t242;
t205 = sin(qJ(4));
t209 = cos(qJ(4));
t142 = t152 * t205 + t155 * t209;
t263 = pkin(1) * t210;
t190 = qJDD(1) * t263;
t197 = qJDD(1) + qJDD(2);
t162 = pkin(2) * t197 - qJD(2) * t242 + t190;
t243 = qJDD(1) * t206;
t266 = pkin(1) * (qJD(2) * t246 + t243);
t149 = t203 * t162 - t202 * t266;
t144 = pkin(3) * t197 + t149;
t150 = t162 * t202 + t203 * t266;
t201 = qJ(1) + qJ(2);
t188 = pkin(9) + qJ(4) + t201;
t182 = sin(t188);
t183 = cos(t188);
t213 = -g(2) * t183 - g(3) * t182 - t142 * qJD(4) + t209 * t144 - t205 * t150;
t193 = qJDD(4) + t197;
t261 = pkin(4) * t193;
t227 = t261 + t213;
t186 = pkin(2) * t203 + pkin(3);
t262 = pkin(2) * t202;
t225 = t186 * t209 - t205 * t262;
t164 = -pkin(4) - t225;
t248 = t205 * t186 + t209 * t262;
t165 = pkin(8) + t248;
t212 = qJD(5) ^ 2;
t253 = t203 * t206;
t223 = pkin(1) * (-t202 * t210 - t253);
t166 = qJD(1) * t223;
t254 = t202 * t206;
t222 = pkin(1) * (t203 * t210 - t254);
t168 = qJD(1) * t222;
t194 = qJD(4) + t198;
t240 = (-t248 * qJD(4) - t166 * t209 + t168 * t205) * t194;
t267 = t164 * t193 + t165 * t212 - t240;
t189 = pkin(2) + t263;
t233 = -pkin(1) * t254 + t203 * t189;
t163 = pkin(3) + t233;
t170 = pkin(1) * t253 + t189 * t202;
t251 = t205 * t163 + t209 * t170;
t255 = t155 * t205;
t214 = g(2) * t182 - (qJD(4) * t152 + t150) * t209 + qJD(4) * t255 - t205 * t144 - g(3) * t183;
t167 = qJD(2) * t223;
t169 = qJD(2) * t222;
t257 = (t251 * qJD(4) - t167 * t209 + t169 * t205) * t194;
t256 = t142 * t194;
t252 = qJDD(3) - g(1);
t250 = -t225 * qJD(4) + t166 * t205 + t168 * t209;
t204 = sin(qJ(5));
t199 = t204 ^ 2;
t208 = cos(qJ(5));
t247 = -t208 ^ 2 + t199;
t245 = qJD(5) * t194;
t244 = qJD(5) * t208;
t195 = sin(t201);
t196 = cos(t201);
t241 = g(2) * t195 - g(3) * t196;
t236 = qJD(1) * (-qJD(2) + t198);
t235 = qJD(2) * (-qJD(1) - t198);
t141 = t152 * t209 - t255;
t139 = -pkin(4) * t194 - t141;
t234 = t139 * t244 - t227 * t204;
t231 = -g(2) * t196 - g(3) * t195;
t177 = qJDD(5) * t204 + t208 * t212;
t178 = qJDD(5) * t208 - t204 * t212;
t230 = 0.2e1 * (t193 * t204 * t208 - t247 * t245) * MDP(12) + (0.2e1 * t194 * t204 * t244 + t193 * t199) * MDP(11) + t177 * MDP(13) + t178 * MDP(14) + t193 * MDP(8);
t229 = t163 * t209 - t170 * t205;
t226 = t190 + t231;
t224 = t197 * MDP(4) + t230;
t220 = pkin(8) * t212 - t256 - t261;
t145 = -pkin(4) - t229;
t146 = pkin(8) + t251;
t219 = t145 * t193 + t146 * t212 + t257;
t218 = -pkin(8) * t193 - t139 * t194 + t214;
t135 = t229 * qJD(4) + t167 * t205 + t169 * t209;
t217 = -qJDD(5) * t146 + (t145 * t194 - t135) * qJD(5);
t216 = -pkin(4) * t245 - pkin(8) * qJDD(5) + qJD(5) * t141;
t215 = -qJDD(5) * t165 + (t164 * t194 + t250) * qJD(5);
t211 = cos(qJ(1));
t207 = sin(qJ(1));
t192 = t194 ^ 2;
t137 = t139 * qJD(5) * t204;
t1 = [qJDD(1) * MDP(1) + (-g(2) * t211 - g(3) * t207) * MDP(2) + (g(2) * t207 - g(3) * t211) * MDP(3) + ((t197 * t210 + t206 * t235) * pkin(1) + t226) * MDP(5) + (((-qJDD(1) - t197) * t206 + t210 * t235) * pkin(1) + t241) * MDP(6) + (t150 * t170 + t155 * t169 + t149 * t233 + t154 * t167 - g(2) * (pkin(1) * t211 + pkin(2) * t196) - g(3) * (pkin(1) * t207 + pkin(2) * t195)) * MDP(7) + (t229 * t193 + t213 - t257) * MDP(9) + (-t135 * t194 - t251 * t193 + t214) * MDP(10) + (t137 + t217 * t204 + (-t219 + t227) * t208) * MDP(16) + (t219 * t204 + t217 * t208 + t234) * MDP(17) + t224; (t236 * t264 + t226) * MDP(5) + ((t210 * t236 - t243) * pkin(1) + t241) * MDP(6) + (-t154 * t166 - t155 * t168 + (t149 * t203 + t150 * t202 + t231) * pkin(2)) * MDP(7) + (t225 * t193 + t213 + t240) * MDP(9) + (-t248 * t193 + t250 * t194 + t214) * MDP(10) + (t137 + t215 * t204 + (t227 - t267) * t208) * MDP(16) + (t267 * t204 + t215 * t208 + t234) * MDP(17) + t224; t178 * MDP(16) - t177 * MDP(17) + t252 * MDP(7); (t213 + t256) * MDP(9) + (t141 * t194 + t214) * MDP(10) + t137 * MDP(16) + t234 * MDP(17) + (t216 * MDP(16) + t220 * MDP(17)) * t204 + ((-t220 + t227) * MDP(16) + t216 * MDP(17)) * t208 + t230; qJDD(5) * MDP(15) + t247 * MDP(12) * t192 + (t193 * MDP(14) + t252 * MDP(16) + t218 * MDP(17)) * t208 + (-t192 * t208 * MDP(11) + t193 * MDP(13) + t218 * MDP(16) - t252 * MDP(17)) * t204;];
tau = t1;
