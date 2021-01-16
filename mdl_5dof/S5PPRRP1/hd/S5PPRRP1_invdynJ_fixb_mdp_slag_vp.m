% Calculate vector of inverse dynamics joint torques for
% S5PPRRP1
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:48:21
% EndTime: 2021-01-15 14:48:24
% DurationCPUTime: 1.31s
% Computational Cost: add. (708->193), mult. (1511->240), div. (0->0), fcn. (1080->10), ass. (0->101)
t185 = sin(pkin(8));
t187 = cos(pkin(8));
t191 = sin(qJ(3));
t193 = cos(qJ(3));
t265 = -t185 * t191 + t187 * t193;
t156 = t265 * qJD(1);
t161 = t185 * t193 + t187 * t191;
t159 = t161 * qJD(3);
t182 = pkin(8) + qJ(3);
t177 = sin(t182);
t186 = sin(pkin(7));
t188 = cos(pkin(7));
t214 = g(1) * t188 + g(2) * t186;
t207 = t214 * t177;
t178 = cos(t182);
t256 = g(3) * t178;
t266 = t207 - t256;
t192 = cos(qJ(4));
t157 = t161 * qJD(1);
t149 = qJD(3) * pkin(6) + t157;
t235 = qJ(5) * qJD(3);
t219 = t149 + t235;
t211 = t219 * t192;
t221 = -g(1) * t186 + g(2) * t188;
t190 = sin(qJ(4));
t183 = t190 ^ 2;
t184 = t192 ^ 2;
t236 = t183 + t184;
t264 = t236 * MDP(15);
t263 = t161 * qJDD(1);
t257 = g(3) * t177;
t200 = t214 * t178 + t257;
t262 = pkin(4) * t183;
t261 = pkin(4) * t192;
t255 = g(3) * t190;
t254 = qJ(5) + pkin(6);
t253 = qJD(3) * pkin(3);
t252 = qJD(4) * pkin(4);
t250 = qJDD(4) * pkin(4);
t249 = t156 * t192;
t246 = t186 * t190;
t245 = t186 * t192;
t242 = t188 * t190;
t241 = t188 * t192;
t195 = qJD(3) ^ 2;
t240 = t192 * t195;
t180 = t192 * qJD(2);
t141 = -t219 * t190 + t180;
t140 = t141 + t252;
t239 = -t141 + t140;
t238 = qJD(4) * t249 + t178 * t255;
t237 = t183 - t184;
t148 = -t156 - t253;
t234 = qJD(3) * t148;
t233 = qJD(3) * t157;
t232 = qJD(3) * t265;
t231 = qJD(3) * t192;
t230 = qJD(4) * t190;
t228 = qJD(3) * qJD(4);
t176 = pkin(3) + t261;
t227 = qJDD(3) * t176;
t226 = qJDD(3) * t190;
t225 = qJDD(3) * t192;
t224 = MDP(11) + MDP(13);
t223 = MDP(12) + MDP(14);
t222 = t190 * t228;
t220 = qJD(4) * t254;
t218 = 0.2e1 * t192 * t228;
t138 = qJDD(3) * pkin(6) + qJD(1) * t232 + t263;
t217 = -qJD(2) * qJD(4) - t138;
t204 = -qJD(1) * t159 + qJDD(1) * t265;
t137 = pkin(4) * t222 + qJDD(5) - t204 - t227;
t216 = t137 - t227;
t215 = t156 * t230 + t157 * t231 + (g(1) * t241 + g(2) * t245) * t177;
t142 = qJD(2) * t190 + t211;
t213 = t140 * t190 - t142 * t192;
t210 = -MDP(5) + t264;
t209 = -qJD(3) * t159 + qJDD(3) * t265;
t194 = qJD(4) ^ 2;
t208 = -0.2e1 * qJDD(3) * pkin(3) + pkin(6) * t194 - t204;
t179 = t192 * qJDD(2);
t206 = -g(1) * (-t178 * t242 + t245) - g(2) * (-t178 * t246 - t241) + t177 * t255 + t179;
t205 = -qJ(5) * qJDD(3) + t217;
t203 = t161 * t194 - t209;
t202 = -pkin(6) * qJDD(4) + (t148 - t253) * qJD(4);
t201 = -0.2e1 * t232 * qJD(4) - qJDD(4) * t161;
t199 = -t207 - t233;
t144 = t149 * t230;
t198 = -g(1) * (-t178 * t241 - t246) - g(2) * (-t178 * t245 + t242) - t190 * qJDD(2) + t144 + t192 * t257;
t197 = qJD(3) * qJD(5) - t205;
t143 = -qJD(3) * t176 + qJD(5) - t156;
t196 = (-qJD(5) - t143) * qJD(3) + t205;
t167 = t254 * t192;
t166 = t254 * t190;
t165 = qJDD(4) * t192 - t190 * t194;
t164 = qJDD(4) * t190 + t192 * t194;
t155 = -qJD(5) * t190 - t192 * t220;
t154 = qJD(5) * t192 - t190 * t220;
t136 = -t144 + (-qJ(5) * t228 + qJDD(2)) * t190 + t197 * t192;
t135 = -qJD(4) * t211 - t197 * t190 + t179 + t250;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (-g(3) + (t185 ^ 2 + t187 ^ 2) * qJDD(1)) * MDP(2) + t209 * MDP(4) + (-t137 * t265 + t143 * t159 - g(3)) * MDP(16) + t224 * (t201 * t190 - t203 * t192) + t223 * (t203 * t190 + t201 * t192) + (-t213 * MDP(16) + t210 * qJD(3)) * t232 + ((-t135 * t190 + t136 * t192 + (-t140 * t192 - t142 * t190) * qJD(4)) * MDP(16) + t210 * qJDD(3)) * t161; (qJDD(2) + t221) * MDP(2) + (-t213 * qJD(4) + t135 * t192 + t136 * t190 + t221) * MDP(16) + t224 * t165 - t223 * t164; qJDD(3) * MDP(3) + (t204 + t233 + t266) * MDP(4) + (-t263 + t200) * MDP(5) + (qJDD(3) * t183 + t190 * t218) * MDP(6) + 0.2e1 * (t190 * t225 - t237 * t228) * MDP(7) + t164 * MDP(8) + t165 * MDP(9) + (t202 * t190 + (-t208 - t256) * t192 + t215) * MDP(11) + (t202 * t192 + (t199 + t208) * t190 + t238) * MDP(12) + (-qJDD(4) * t166 + (-t216 - t256) * t192 + (t155 + (t143 + (-t176 - t261) * qJD(3)) * t190) * qJD(4) + t215) * MDP(13) + (-qJDD(4) * t167 + (t143 * t192 - t154 + (-t176 * t192 + t262) * qJD(3)) * qJD(4) + (t199 + t216) * t190 + t238) * MDP(14) + ((-qJD(4) * t140 + qJDD(3) * t167 + t136) * t192 + (-qJD(4) * t142 + qJDD(3) * t166 - t135) * t190 + (t154 * t192 - t155 * t190 - t236 * t156 + (t166 * t192 - t167 * t190) * qJD(4)) * qJD(3) - t200) * MDP(15) + (t136 * t167 - t135 * t166 - t137 * t176 - g(3) * (t176 * t178 + t177 * t254) + (pkin(4) * t230 - t157) * t143 + (t154 - t249) * t142 + (t156 * t190 + t155) * t140 + t214 * (t176 * t177 - t178 * t254)) * MDP(16); -t190 * MDP(6) * t240 + t237 * t195 * MDP(7) + MDP(8) * t226 + MDP(9) * t225 + qJDD(4) * MDP(10) + ((-t138 - t234) * t190 + t206) * MDP(11) + ((-t149 * t190 + t180) * qJD(4) + (t217 - t234) * t192 + t198) * MDP(12) + (0.2e1 * t250 + (t142 - t211) * qJD(4) + (pkin(4) * t240 + t196) * t190 + t206) * MDP(13) + (-t195 * t262 + (t190 * t235 + t141) * qJD(4) + t196 * t192 + t198) * MDP(14) + (-pkin(4) * t226 + (t239 - t252) * t231) * MDP(15) + (t239 * t142 + (t135 + t221 * t192 + (-t143 * qJD(3) + t200) * t190) * pkin(4)) * MDP(16); (0.2e1 * t222 - t225) * MDP(13) + (t218 + t226) * MDP(14) + (t213 * qJD(3) + t137 - t266) * MDP(16) - t195 * t264;];
tau = t1;
