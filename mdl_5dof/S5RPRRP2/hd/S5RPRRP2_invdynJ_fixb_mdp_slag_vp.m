% Calculate vector of inverse dynamics joint torques for
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:35
% EndTime: 2020-01-03 11:45:38
% DurationCPUTime: 1.25s
% Computational Cost: add. (1031->182), mult. (1732->232), div. (0->0), fcn. (947->12), ass. (0->101)
t201 = qJ(1) + pkin(8);
t194 = qJ(3) + t201;
t188 = sin(t194);
t189 = cos(t194);
t229 = g(2) * t189 + g(3) * t188;
t208 = sin(qJ(3));
t211 = cos(qJ(3));
t205 = cos(pkin(8));
t190 = pkin(1) * t205 + pkin(2);
t176 = t190 * qJD(1);
t204 = sin(pkin(8));
t268 = pkin(1) * t204;
t272 = qJD(3) * t176 + qJDD(1) * t268;
t243 = qJD(1) * t268;
t273 = -qJD(3) * t243 + t190 * qJDD(1);
t230 = -t272 * t208 + t273 * t211;
t218 = t229 - t230;
t206 = -qJ(5) - pkin(7);
t199 = qJDD(1) + qJDD(3);
t262 = t199 * pkin(3);
t223 = t262 - t218;
t271 = g(2) * t188 - g(3) * t189;
t234 = t190 * t211 - t208 * t268;
t248 = t208 * t190 + t211 * t268;
t207 = sin(qJ(4));
t202 = t207 ^ 2;
t210 = cos(qJ(4));
t203 = t210 ^ 2;
t270 = -t202 - t203;
t269 = -t273 * t208 - t272 * t211;
t200 = qJD(1) + qJD(3);
t267 = pkin(3) * t200;
t266 = pkin(4) * t210;
t265 = g(1) * t210;
t261 = qJD(4) * pkin(4);
t260 = pkin(7) * qJDD(4);
t154 = t211 * t176 - t208 * t243;
t259 = t154 * t200;
t155 = t176 * t208 + t211 * t243;
t258 = t155 * t200;
t158 = t248 * qJD(3);
t257 = t158 * t200;
t255 = t199 * t207;
t254 = t207 * t210;
t160 = pkin(7) + t248;
t253 = -qJ(5) - t160;
t252 = qJDD(2) - g(1);
t237 = -t206 * t200 + t155;
t144 = t210 * qJD(2) - t237 * t207;
t143 = t144 + t261;
t251 = t143 - t144;
t191 = pkin(3) + t266;
t250 = t188 * t191 + t189 * t206;
t247 = t202 - t203;
t246 = MDP(16) * t154;
t244 = qJD(4) * t207;
t141 = pkin(7) * t199 - t269;
t216 = qJ(5) * t199 + qJD(2) * qJD(4) + qJD(5) * t200 + t141;
t224 = t237 * qJD(4);
t136 = (qJDD(2) - t224) * t207 + t216 * t210;
t241 = t136 * t210 - t271;
t240 = t200 * t244;
t239 = t154 - t267;
t238 = qJD(4) * t206;
t236 = -t188 * t206 + t189 * t191;
t233 = qJD(4) * t253;
t150 = -t154 - t267;
t231 = t150 * qJD(4) * t210 - t223 * t207;
t159 = -pkin(3) - t234;
t209 = sin(qJ(1));
t212 = cos(qJ(1));
t227 = -g(2) * t212 - g(3) * t209;
t213 = qJD(4) ^ 2;
t169 = qJDD(4) * t207 + t210 * t213;
t170 = qJDD(4) * t210 - t207 * t213;
t226 = 0.2e1 * (-t247 * t200 * qJD(4) + t199 * t254) * MDP(9) + (t199 * t202 + 0.2e1 * t210 * t240) * MDP(8) + t169 * MDP(10) + t170 * MDP(11) + t199 * MDP(5);
t145 = qJD(2) * t207 + t237 * t210;
t225 = t143 * t207 - t145 * t210;
t221 = pkin(7) * t213 - t258 - t262;
t220 = t159 * t199 + t160 * t213 + t257;
t219 = -t150 * t200 - t141 + t271;
t157 = t234 * qJD(3);
t217 = -qJDD(4) * t160 + (t159 * t200 - t157) * qJD(4);
t137 = pkin(4) * t240 - t191 * t199 + qJDD(5) - t230;
t215 = t271 + t269;
t198 = t200 ^ 2;
t197 = t210 * qJ(5);
t195 = t210 * qJD(5);
t193 = t210 * qJDD(2);
t179 = pkin(7) * t210 + t197;
t178 = t206 * t207;
t162 = -qJD(5) * t207 + t210 * t238;
t161 = t207 * t238 + t195;
t153 = t160 * t210 + t197;
t152 = t253 * t207;
t147 = t150 * t244;
t146 = -t191 * t200 + qJD(5) - t154;
t139 = (-qJD(5) - t157) * t207 + t210 * t233;
t138 = t157 * t210 + t207 * t233 + t195;
t135 = qJDD(4) * pkin(4) - t216 * t207 - t210 * t224 + t193;
t1 = [qJDD(1) * MDP(1) + t227 * MDP(2) + (g(2) * t209 - g(3) * t212) * MDP(3) + (t227 + (t204 ^ 2 + t205 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t234 * t199 - t218 - t257) * MDP(6) + (-t157 * t200 - t248 * t199 + t215) * MDP(7) + (t147 + t217 * t207 + (-t220 + t223) * t210) * MDP(13) + (t220 * t207 + t217 * t210 + t231) * MDP(14) + ((t138 * t200 + t153 * t199 + (-t152 * t200 - t143) * qJD(4)) * t210 + (-t139 * t200 - t152 * t199 - t135 + (-t153 * t200 - t145) * qJD(4)) * t207 + t241) * MDP(15) + (t136 * t153 + t145 * t138 + t135 * t152 + t143 * t139 + t137 * (t159 - t266) + t146 * (pkin(4) * t244 + t158) - g(2) * (pkin(2) * cos(t201) + t212 * pkin(1) + t236) - g(3) * (pkin(2) * sin(t201) + t209 * pkin(1) + t250)) * MDP(16) + t226; t252 * MDP(4) + t170 * MDP(13) - t169 * MDP(14) + (-t225 * qJD(4) + t135 * t210 + t136 * t207 - g(1)) * MDP(16); (-t218 + t258) * MDP(6) + (t215 + t259) * MDP(7) + t147 * MDP(13) + t231 * MDP(14) + (t270 * t259 + t241) * MDP(15) + (-g(2) * t236 - g(3) * t250 + t135 * t178 + t136 * t179 - t137 * t191 + t143 * t162 + t145 * t161 - t146 * t155) * MDP(16) + (-MDP(13) * t260 + t221 * MDP(14) + (-t162 * t200 - t178 * t199 - t135) * MDP(15) + t143 * t246 + (t239 * MDP(13) + (-t179 * t200 - t145) * MDP(15) + pkin(4) * t146 * MDP(16)) * qJD(4)) * t207 + ((-t221 + t223) * MDP(13) - MDP(14) * t260 + (t161 * t200 + t179 * t199) * MDP(15) - t145 * t246 + (t239 * MDP(14) + (-t178 * t200 - t143) * MDP(15)) * qJD(4)) * t210 + t226; MDP(10) * t255 + t210 * t199 * MDP(11) + qJDD(4) * MDP(12) + (t219 * t207 + t193 - t265) * MDP(13) + (-t252 * t207 + t219 * t210) * MDP(14) + (-pkin(4) * t255 + (t251 - t261) * t210 * t200) * MDP(15) + (t251 * t145 + (-t265 + t135 + (-t146 * t200 + t271) * t207) * pkin(4)) * MDP(16) + (-MDP(8) * t254 + t247 * MDP(9)) * t198; (t225 * t200 + t137 + t229) * MDP(16) + t270 * MDP(15) * t198;];
tau = t1;
