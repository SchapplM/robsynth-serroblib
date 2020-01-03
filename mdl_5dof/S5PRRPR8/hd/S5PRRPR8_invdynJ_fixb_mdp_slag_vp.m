% Calculate vector of inverse dynamics joint torques for
% S5PRRPR8
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
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PRRPR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:46
% EndTime: 2019-12-31 17:42:48
% DurationCPUTime: 0.92s
% Computational Cost: add. (850->162), mult. (1559->229), div. (0->0), fcn. (1159->14), ass. (0->94)
t207 = qJ(2) + qJ(3);
t199 = pkin(9) + t207;
t191 = sin(t199);
t209 = sin(pkin(8));
t211 = cos(pkin(8));
t234 = g(1) * t211 + g(2) * t209;
t262 = t234 * t191;
t200 = sin(t207);
t201 = cos(t207);
t265 = -g(3) * t201 + t234 * t200;
t217 = cos(qJ(2));
t198 = t217 * qJDD(1);
t214 = sin(qJ(2));
t242 = qJD(1) * qJD(2);
t171 = qJDD(2) * pkin(2) - t214 * t242 + t198;
t216 = cos(qJ(3));
t166 = t216 * t171;
t203 = qJDD(2) + qJDD(3);
t213 = sin(qJ(3));
t189 = qJD(2) * pkin(2) + qJD(1) * t217;
t228 = -qJDD(1) * t214 - t217 * t242;
t223 = qJD(3) * t189 - t228;
t245 = qJD(1) * t214;
t237 = qJD(3) * t245;
t140 = pkin(3) * t203 - t223 * t213 - t216 * t237 + t166;
t188 = t213 * t237;
t146 = t171 * t213 + t223 * t216 - t188;
t208 = sin(pkin(9));
t210 = cos(pkin(9));
t135 = t140 * t210 - t146 * t208;
t192 = cos(t199);
t263 = -pkin(4) * t203 + g(3) * t192 - t135;
t204 = qJD(2) + qJD(3);
t196 = pkin(2) * t216 + pkin(3);
t253 = t208 * t213;
t230 = -pkin(2) * t253 + t196 * t210;
t164 = -pkin(4) - t230;
t251 = t210 * t213;
t247 = pkin(2) * t251 + t208 * t196;
t165 = pkin(7) + t247;
t218 = qJD(5) ^ 2;
t175 = t213 * t217 + t214 * t216;
t169 = t175 * qJD(1);
t174 = -t213 * t214 + t216 * t217;
t170 = t174 * qJD(1);
t255 = pkin(2) * qJD(3);
t249 = t169 * t210 + t170 * t208 - (t208 * t216 + t251) * t255;
t261 = -t164 * t203 - t165 * t218 + t249 * t204;
t260 = pkin(2) * t203;
t163 = t189 * t213 + t216 * t245;
t254 = t163 * t208;
t252 = t210 * t163;
t250 = qJDD(1) - g(3);
t136 = t208 * t140 + t210 * t146;
t248 = t169 * t208 - t170 * t210 + (t210 * t216 - t253) * t255;
t212 = sin(qJ(5));
t205 = t212 ^ 2;
t215 = cos(qJ(5));
t246 = -t215 ^ 2 + t205;
t244 = qJD(5) * t204;
t243 = t215 * qJD(5);
t162 = t216 * t189 - t213 * t245;
t159 = pkin(3) * t204 + t162;
t147 = t159 * t210 - t254;
t144 = -pkin(4) * t204 - t147;
t239 = t144 * t243 + t263 * t212;
t238 = t144 * qJD(5) * t212 + t262 * t215;
t235 = g(3) * t200 + t234 * t201 + t188;
t180 = qJDD(5) * t212 + t215 * t218;
t181 = qJDD(5) * t215 - t212 * t218;
t233 = 0.2e1 * (t203 * t212 * t215 - t246 * t244) * MDP(10) + (0.2e1 * t204 * t212 * t243 + t203 * t205) * MDP(9) + t180 * MDP(11) + t181 * MDP(12) + t203 * MDP(5);
t231 = t166 + t265;
t229 = g(1) * t209 - g(2) * t211 - qJDD(4);
t156 = t204 * t174;
t157 = t204 * t175;
t137 = t156 * t208 + t157 * t210;
t153 = -t174 * t210 + t175 * t208;
t154 = t174 * t208 + t175 * t210;
t227 = t137 * t204 + t153 * t203 + t154 * t218;
t149 = t162 * t208 + t252;
t193 = pkin(3) * t208 + pkin(7);
t194 = -pkin(3) * t210 - pkin(4);
t226 = t149 * t204 - t193 * t218 - t194 * t203;
t138 = t156 * t210 - t157 * t208;
t225 = -qJD(5) * t138 - qJDD(5) * t154 + t153 * t244;
t150 = t162 * t210 - t254;
t224 = qJD(5) * t150 - qJDD(5) * t193 + t194 * t244;
t222 = -qJDD(5) * t165 + (t164 * t204 - t248) * qJD(5);
t221 = (-pkin(2) * t204 - t189) * qJD(3) + t228;
t220 = -pkin(7) * t203 + g(3) * t191 - t144 * t204 + t234 * t192 - t136;
t219 = qJD(2) ^ 2;
t202 = t204 ^ 2;
t148 = t208 * t159 + t252;
t1 = [t250 * MDP(1) + (qJDD(2) * t217 - t214 * t219) * MDP(3) + (-qJDD(2) * t214 - t217 * t219) * MDP(4) + (-t157 * t204 + t174 * t203) * MDP(6) + (-t156 * t204 - t175 * t203) * MDP(7) + (-t135 * t153 + t136 * t154 - t137 * t147 + t138 * t148 - g(3)) * MDP(8) + (-t227 * MDP(14) + t225 * MDP(15)) * t215 + (t225 * MDP(14) + t227 * MDP(15)) * t212; qJDD(2) * MDP(2) + (-g(3) * t217 + t234 * t214 + t198) * MDP(3) + (-t250 * t214 + t234 * t217) * MDP(4) + (t169 * t204 + (-t237 + t260) * t216 + t221 * t213 + t231) * MDP(6) + (t170 * t204 + (-t171 - t260) * t213 + t221 * t216 + t235) * MDP(7) + (t136 * t247 + t135 * t230 - g(3) * (pkin(2) * t217 + pkin(3) * t201) - t234 * (-pkin(2) * t214 - pkin(3) * t200) + t248 * t148 + t249 * t147) * MDP(8) + (t222 * t212 + (-t263 + t261) * t215 + t238) * MDP(14) + (t222 * t215 + (-t261 - t262) * t212 + t239) * MDP(15) + t233; (t163 * t204 + t231) * MDP(6) + (t162 * t204 + t235) * MDP(7) + t238 * MDP(14) + t239 * MDP(15) + (-MDP(6) * t237 - t223 * MDP(7)) * t216 + (-t223 * MDP(6) - t171 * MDP(7)) * t213 + ((t226 - t263) * MDP(14) + t224 * MDP(15)) * t215 + (t224 * MDP(14) + (-t226 - t262) * MDP(15)) * t212 + t233 + (t147 * t149 - t148 * t150 + (t135 * t210 + t136 * t208 + t265) * pkin(3)) * MDP(8); t181 * MDP(14) - t180 * MDP(15) - t229 * MDP(8); qJDD(5) * MDP(13) + t246 * MDP(10) * t202 + (t203 * MDP(12) - t229 * MDP(14) + t220 * MDP(15)) * t215 + (-t202 * t215 * MDP(9) + t203 * MDP(11) + t220 * MDP(14) + t229 * MDP(15)) * t212;];
tau = t1;
