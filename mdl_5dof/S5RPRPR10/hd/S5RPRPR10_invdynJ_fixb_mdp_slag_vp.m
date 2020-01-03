% Calculate vector of inverse dynamics joint torques for
% S5RPRPR10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPR10_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:08
% EndTime: 2019-12-31 18:26:10
% DurationCPUTime: 0.87s
% Computational Cost: add. (932->165), mult. (1380->220), div. (0->0), fcn. (758->10), ass. (0->90)
t231 = qJD(1) - qJD(3);
t209 = -pkin(1) - pkin(2);
t237 = qJ(2) * qJD(1);
t256 = -qJD(3) * t237 + t209 * qJDD(1) + qJDD(2);
t180 = t209 * qJD(1) + qJD(2);
t233 = (qJD(1) * qJD(2));
t223 = qJD(3) * t180 + t233;
t234 = (qJ(2) * qJDD(1));
t255 = -t223 - t234;
t201 = sin(pkin(8));
t202 = cos(pkin(8));
t204 = sin(qJ(3));
t207 = cos(qJ(3));
t167 = t201 * t207 + t202 * t204;
t254 = t231 * t167;
t166 = t201 * t204 - t202 * t207;
t253 = t231 * t166;
t189 = t207 * t209;
t252 = -qJ(2) * t204 + t189;
t208 = cos(qJ(1));
t205 = sin(qJ(1));
t243 = t205 * t204;
t169 = -t208 * t207 - t243;
t244 = t204 * t208;
t170 = -t205 * t207 + t244;
t251 = g(1) * t170 - g(2) * t169;
t248 = (pkin(1) * qJDD(1));
t157 = t204 * t180 + t207 * t237;
t246 = t157 * t201;
t245 = t202 * t157;
t242 = qJDD(4) + g(3);
t197 = qJDD(1) - qJDD(3);
t221 = t256 * t207;
t212 = t255 * t204 + t221;
t138 = -pkin(3) * t197 + t212;
t144 = t256 * t204 - t255 * t207;
t134 = t201 * t138 + t202 * t144;
t171 = -pkin(3) + t252;
t175 = qJ(2) * t207 + t204 * t209;
t241 = t201 * t171 + t202 * t175;
t240 = t208 * pkin(1) + t205 * qJ(2);
t239 = g(1) * t205 - g(2) * t208;
t203 = sin(qJ(5));
t199 = t203 ^ 2;
t206 = cos(qJ(5));
t238 = -t206 ^ 2 + t199;
t235 = qJD(5) * t206;
t232 = qJ(3) + pkin(8);
t230 = 2 * t233;
t228 = cos(t232);
t227 = qJDD(2) - t248;
t226 = g(1) * t208 + g(2) * t205;
t156 = t207 * t180 - t204 * t237;
t133 = t138 * t202 - t144 * t201;
t151 = -pkin(3) * t231 + t156;
t139 = t151 * t202 - t246;
t210 = qJD(5) ^ 2;
t225 = t166 * t197 - t167 * t210;
t224 = t171 * t202 - t175 * t201;
t190 = sin(t232);
t158 = t190 * t208 - t205 * t228;
t159 = t205 * t190 + t208 * t228;
t222 = -pkin(4) * t197 + g(1) * t158 + g(2) * t159 + t133;
t135 = pkin(4) * t231 - t139;
t220 = pkin(7) * t197 + g(1) * t159 - g(2) * t158 + t135 * t231 - t134;
t219 = -t253 * qJD(5) - qJDD(5) * t167;
t154 = qJD(2) * t207 + t252 * qJD(3);
t155 = -t204 * qJD(2) - t175 * qJD(3);
t143 = t154 * t202 + t155 * t201;
t147 = pkin(4) - t224;
t148 = -pkin(7) + t241;
t218 = -qJDD(5) * t148 + (-t147 * t231 - t135 - t143) * qJD(5);
t146 = t156 * t202 - t246;
t181 = pkin(3) * t201 + pkin(7);
t182 = -pkin(3) * t202 - pkin(4);
t217 = -qJDD(5) * t181 + (-t182 * t231 + t135 + t146) * qJD(5);
t142 = t154 * t201 - t155 * t202;
t216 = -t142 * t231 - t147 * t197 + t148 * t210 + t222;
t145 = t156 * t201 + t245;
t215 = -t145 * t231 - t181 * t210 + t182 * t197 + t222;
t176 = qJDD(5) * t203 + t206 * t210;
t177 = qJDD(5) * t206 - t203 * t210;
t214 = (-0.2e1 * t203 * t231 * t235 - t197 * t199) * MDP(11) + 0.2e1 * (qJD(5) * t231 * t238 - t197 * t203 * t206) * MDP(12) + t176 * MDP(13) + t177 * MDP(14) - t197 * MDP(7);
t213 = -g(1) * t169 - g(2) * t170 - t144;
t211 = qJD(1) ^ 2;
t196 = t231 ^ 2;
t192 = t208 * qJ(2);
t188 = pkin(3) * t207 + pkin(2);
t140 = t201 * t151 + t245;
t1 = [qJDD(1) * MDP(1) + t239 * MDP(2) + t226 * MDP(3) + (-qJDD(2) + t239 + (2 * t248)) * MDP(4) + (-t226 + t230 + (2 * t234)) * MDP(5) + (-t227 * pkin(1) - g(1) * (-pkin(1) * t205 + t192) - g(2) * t240 + (t230 + t234) * qJ(2)) * MDP(6) + (-t155 * t231 - t189 * t197 + ((qJDD(1) + t197) * qJ(2) + t223) * t204 - t221 - t251) * MDP(8) + (t154 * t231 + t175 * t197 - t213) * MDP(9) + (t134 * t241 + t140 * t143 + t133 * t224 - t139 * t142 - g(1) * (pkin(3) * t244 + t192 + (-pkin(1) - t188) * t205) - g(2) * (pkin(3) * t243 + t188 * t208 + t240)) * MDP(10) + (t218 * t203 - t216 * t206) * MDP(16) + (t216 * t203 + t218 * t206) * MDP(17) - t214; -qJDD(1) * MDP(4) - t211 * MDP(5) + (-qJ(2) * t211 + t227 - t239) * MDP(6) + (-t133 * t166 + t134 * t167 + t254 * t139 + t253 * t140 - t239) * MDP(10) + (-MDP(8) * t207 + MDP(9) * t204) * t197 + (t225 * MDP(16) + t219 * MDP(17)) * t206 + (t219 * MDP(16) - t225 * MDP(17)) * t203 - ((qJD(5) * t166 * t203 + t254 * t206) * MDP(16) + (t166 * t235 - t254 * t203) * MDP(17) + (MDP(8) * t204 + MDP(9) * t207) * t231) * t231; (-t157 * t231 + t212 + t251) * MDP(8) + (-t156 * t231 + t213) * MDP(9) + (t139 * t145 - t140 * t146 + (t133 * t202 + t134 * t201 + t251) * pkin(3)) * MDP(10) + (t217 * t203 + t215 * t206) * MDP(16) + (-t215 * t203 + t217 * t206) * MDP(17) + t214; t242 * MDP(10) + t177 * MDP(16) - t176 * MDP(17); qJDD(5) * MDP(15) + t238 * MDP(12) * t196 + (-t197 * MDP(14) + t242 * MDP(16) + t220 * MDP(17)) * t206 + (-t196 * t206 * MDP(11) - t197 * MDP(13) + t220 * MDP(16) - t242 * MDP(17)) * t203;];
tau = t1;
