% Calculate vector of inverse dynamics joint torques for
% S5PPRRP2
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
%   see S5PPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:09:15
% EndTime: 2019-12-05 15:09:18
% DurationCPUTime: 1.22s
% Computational Cost: add. (721->174), mult. (1497->214), div. (0->0), fcn. (1089->10), ass. (0->100)
t189 = sin(pkin(8));
t191 = cos(pkin(8));
t194 = sin(qJ(3));
t196 = cos(qJ(3));
t271 = -t189 * t194 + t191 * t196;
t160 = t271 * qJD(1);
t165 = t189 * t196 + t191 * t194;
t163 = t165 * qJD(3);
t236 = qJD(3) * t271;
t254 = qJDD(3) * pkin(6);
t265 = t165 * qJDD(1);
t272 = qJD(1) * t236 + qJD(2) * qJD(4) + t254 + t265;
t186 = pkin(8) + qJ(3);
t182 = cos(t186);
t258 = g(3) * t182;
t181 = sin(t186);
t190 = sin(pkin(7));
t192 = cos(pkin(7));
t214 = g(1) * t192 + g(2) * t190;
t267 = t214 * t181;
t203 = -t258 + t267;
t193 = sin(qJ(4));
t195 = cos(qJ(4));
t209 = pkin(4) * t195 + qJ(5) * t193 + pkin(3);
t270 = t209 * qJD(3);
t161 = t165 * qJD(1);
t154 = qJD(3) * pkin(6) + t161;
t252 = t154 * t193;
t148 = qJD(2) * t195 - t252;
t269 = qJD(5) - t148;
t268 = t209 * qJDD(3);
t145 = -qJD(4) * pkin(4) + t269;
t251 = t154 * t195;
t149 = qJD(2) * t193 + t251;
t146 = qJD(4) * qJ(5) + t149;
t208 = -qJD(1) * t163 + t271 * qJDD(1);
t213 = pkin(4) * t193 - qJ(5) * t195;
t159 = t213 * qJD(4) - qJD(5) * t193;
t238 = qJD(3) * t159;
t140 = -t208 + t238 - t268;
t266 = -t140 + t268;
t253 = qJDD(4) * pkin(4);
t264 = qJDD(5) - t253;
t259 = g(3) * t181;
t263 = t214 * t182 + t259;
t197 = qJD(4) ^ 2;
t262 = pkin(6) * t197;
t257 = qJD(3) * pkin(3);
t256 = pkin(6) * qJDD(4);
t248 = t190 * t193;
t247 = t190 * t195;
t244 = t192 * t193;
t243 = t192 * t195;
t242 = t193 * t195;
t241 = t159 - t161;
t187 = t193 ^ 2;
t188 = t195 ^ 2;
t240 = t187 - t188;
t239 = t187 + t188;
t237 = qJD(3) * t161;
t235 = qJD(3) * t193;
t231 = qJD(3) * qJD(4);
t229 = qJDD(3) * t193;
t228 = qJDD(3) * t195;
t227 = qJDD(4) * qJ(5);
t226 = MDP(11) + MDP(13);
t225 = MDP(12) - MDP(15);
t198 = qJD(3) ^ 2;
t224 = t198 * t242;
t223 = t193 * qJDD(2) + t272 * t195;
t222 = -g(1) * t190 + g(2) * t192;
t153 = -t160 - t257;
t221 = t153 - t257;
t219 = t148 + t252;
t147 = -t160 - t270;
t218 = t147 - t270;
t217 = qJD(4) * t251 - t195 * qJDD(2) + t272 * t193;
t216 = t193 * t160 * qJD(4) + t195 * t237 + (g(1) * t243 + g(2) * t247) * t181;
t215 = t258 + t262;
t212 = t145 * t193 + t146 * t195;
t211 = t239 * MDP(14) - MDP(5);
t210 = -qJD(3) * t163 + qJDD(3) * t271;
t156 = t182 * t247 - t244;
t158 = t182 * t243 + t248;
t207 = g(1) * t158 + g(2) * t156 - t223;
t206 = t165 * t197 - t210;
t205 = 0.2e1 * qJDD(3) * pkin(3) + t208 - t215;
t204 = -0.2e1 * t236 * qJD(4) - qJDD(4) * t165;
t155 = t182 * t248 + t243;
t157 = t182 * t244 - t247;
t202 = g(1) * t157 + g(2) * t155 + t193 * t259 - t217;
t201 = qJD(4) * t149 + t202;
t138 = t227 + (qJD(5) - t252) * qJD(4) + t223;
t139 = t217 + t264;
t200 = t138 * t195 + t139 * t193 + (t145 * t195 - t146 * t193) * qJD(4);
t199 = t200 - t263;
t171 = qJDD(4) * t195 - t193 * t197;
t170 = qJDD(4) * t193 + t195 * t197;
t166 = t213 * qJD(3);
t1 = [(qJDD(1) - g(3)) * MDP(1) + (-g(3) + (t189 ^ 2 + t191 ^ 2) * qJDD(1)) * MDP(2) + t210 * MDP(4) + (-t140 * t271 + t147 * t163 - g(3)) * MDP(16) + t226 * (t204 * t193 - t206 * t195) + t225 * (t206 * t193 + t204 * t195) + (t212 * MDP(16) + t211 * qJD(3)) * t236 + (t200 * MDP(16) + t211 * qJDD(3)) * t165; (qJDD(2) + t222) * MDP(2) + (t212 * qJD(4) + t138 * t193 - t139 * t195 + t222) * MDP(16) + t226 * t171 - t225 * t170; qJDD(3) * MDP(3) + (t203 + t208 + t237) * MDP(4) + (-t265 + t263) * MDP(5) + (qJDD(3) * t187 + 0.2e1 * t231 * t242) * MDP(6) + 0.2e1 * (t193 * t228 - t240 * t231) * MDP(7) + t170 * MDP(8) + t171 * MDP(9) + ((t221 * qJD(4) - t256) * t193 + t205 * t195 + t216) * MDP(11) + ((-t256 + (t160 + t221) * qJD(4)) * t195 + (-t205 - t237 - t267) * t193) * MDP(12) + ((t218 * qJD(4) - t256) * t193 + (-t215 - t238 + t266) * t195 + t216) * MDP(13) + (t199 + (-qJD(3) * t160 + t254) * t239) * MDP(14) + ((t256 + (-t160 - t218) * qJD(4)) * t195 + (-t241 * qJD(3) + t203 - t262 + t266) * t193) * MDP(15) + (t199 * pkin(6) + t241 * t147 - t212 * t160 + (-t140 + t203) * t209) * MDP(16); -MDP(6) * t224 + t240 * t198 * MDP(7) + MDP(8) * t229 + MDP(9) * t228 + qJDD(4) * MDP(10) + (-t153 * t235 + t201) * MDP(11) + ((-qJD(3) * t153 + t259) * t195 + t219 * qJD(4) + t207) * MDP(12) + (0.2e1 * t253 - qJDD(5) + (-t147 * t193 + t166 * t195) * qJD(3) + t201) * MDP(13) - t213 * qJDD(3) * MDP(14) + (-t195 * t259 + 0.2e1 * t227 + (t147 * t195 + t166 * t193) * qJD(3) + (0.2e1 * qJD(5) - t219) * qJD(4) - t207) * MDP(15) + (t138 * qJ(5) - t139 * pkin(4) - t147 * t166 - t145 * t149 - g(1) * (-pkin(4) * t157 + qJ(5) * t158) - g(2) * (-pkin(4) * t155 + qJ(5) * t156) + t213 * t259 + t269 * t146) * MDP(16); (-qJDD(4) - t224) * MDP(13) + MDP(14) * t229 + (-t187 * t198 - t197) * MDP(15) + (-qJD(4) * t146 + t147 * t235 - t202 + t264) * MDP(16);];
tau = t1;
