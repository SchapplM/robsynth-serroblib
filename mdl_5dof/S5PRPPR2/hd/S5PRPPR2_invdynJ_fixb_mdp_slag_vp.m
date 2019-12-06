% Calculate vector of inverse dynamics joint torques for
% S5PRPPR2
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PRPPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:47
% EndTime: 2019-12-05 15:24:51
% DurationCPUTime: 1.56s
% Computational Cost: add. (637->184), mult. (1302->253), div. (0->0), fcn. (1010->14), ass. (0->102)
t205 = sin(pkin(9));
t208 = cos(pkin(9));
t251 = t208 * MDP(6);
t271 = -MDP(7) * t205 + t251;
t212 = sin(qJ(2));
t214 = cos(qJ(2));
t242 = qJD(1) * qJD(2);
t269 = qJDD(1) * t212 + t214 * t242;
t204 = qJ(2) + pkin(8);
t196 = sin(t204);
t198 = cos(t204);
t207 = sin(pkin(7));
t210 = cos(pkin(7));
t228 = g(1) * t210 + g(2) * t207;
t218 = -g(3) * t198 + t196 * t228;
t247 = qJD(1) * t214;
t179 = qJD(2) * pkin(2) + t247;
t209 = cos(pkin(8));
t248 = qJD(1) * t212;
t186 = t209 * t248;
t206 = sin(pkin(8));
t156 = t206 * t179 + t186;
t154 = qJD(2) * qJ(4) + t156;
t149 = t208 * qJD(3) - t154 * t205;
t150 = t205 * qJD(3) + t208 * t154;
t267 = -t149 * t205 + t150 * t208;
t249 = t205 ^ 2 + t208 ^ 2;
t266 = qJD(2) * t249;
t211 = sin(qJ(5));
t213 = cos(qJ(5));
t173 = t205 * t213 + t208 * t211;
t166 = t173 * qJD(5);
t185 = t206 * t248;
t163 = t209 * t247 - t185;
t243 = qJD(4) - t163;
t264 = qJD(5) ^ 2;
t263 = pkin(2) * t209;
t260 = g(3) * t196;
t188 = pkin(2) * t206 + qJ(4);
t258 = pkin(6) + t188;
t256 = pkin(6) * qJDD(2);
t253 = t198 * t207;
t252 = t198 * t210;
t250 = qJDD(1) - g(3);
t199 = t214 * qJDD(1);
t169 = qJDD(2) * pkin(2) - t212 * t242 + t199;
t144 = t206 * t169 + t269 * t209;
t139 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t144;
t136 = t205 * qJDD(3) + t208 * t139;
t245 = qJD(2) * t211;
t244 = qJD(2) * t213;
t239 = qJDD(2) * t211;
t238 = qJDD(2) * t213;
t171 = t205 * t211 - t213 * t208;
t237 = qJDD(5) * t171;
t236 = qJDD(5) * t173;
t232 = t208 * t244;
t235 = qJD(5) * t232 + t205 * t238 + t208 * t239;
t234 = -pkin(4) * t208 - pkin(3);
t233 = t205 * t245;
t230 = -g(1) * t207 + g(2) * t210;
t155 = t179 * t209 - t185;
t229 = qJDD(2) * t249;
t143 = t169 * t209 - t269 * t206;
t227 = qJD(4) - t155;
t184 = t208 * t238;
t226 = -t205 * t239 + t184;
t192 = t208 * qJDD(3);
t135 = -t139 * t205 + t192;
t225 = -t135 * t205 + t136 * t208;
t167 = t258 * t205;
t168 = t258 * t208;
t224 = -t167 * t213 - t168 * t211;
t223 = -t167 * t211 + t168 * t213;
t172 = t206 * t214 + t209 * t212;
t170 = t206 * t212 - t209 * t214;
t221 = qJDD(4) - t143;
t165 = t171 * qJD(5);
t219 = -t205 * t244 - t208 * t245;
t216 = -g(3) * t214 + t212 * t228;
t215 = qJD(2) ^ 2;
t203 = pkin(9) + qJ(5);
t197 = cos(t203);
t195 = sin(t203);
t190 = -pkin(3) - t263;
t175 = t234 - t263;
t164 = t173 * qJD(2);
t162 = t170 * qJD(2);
t160 = -t232 + t233;
t159 = t206 * t247 + t186;
t158 = t172 * qJD(2);
t153 = -qJD(2) * pkin(3) + t227;
t151 = qJD(2) * t234 + t227;
t148 = -qJD(5) * t166 - t237;
t147 = -qJD(5) * t165 + t236;
t146 = qJD(2) * t166 - t226;
t145 = -qJD(5) * t233 + t235;
t140 = -qJDD(2) * pkin(3) + t221;
t137 = qJDD(2) * t234 + t221;
t134 = t208 * t256 + t136;
t133 = t192 + (-t139 - t256) * t205;
t1 = [t250 * MDP(1) + (qJDD(2) * t214 - t212 * t215) * MDP(3) + (-qJDD(2) * t212 - t214 * t215) * MDP(4) + (-t143 * t170 + t144 * t172 - t155 * t158 - t156 * t162 - g(3)) * MDP(5) + (-t162 * t266 + t172 * t229) * MDP(8) + (t140 * t170 + t153 * t158 - t267 * t162 + t225 * t172 - g(3)) * MDP(9) + (t170 * t146 + t158 * t160 + t162 * t166 + (t171 * t264 - t236) * t172) * MDP(15) + (t170 * t145 + t158 * t164 - t162 * t165 + (t173 * t264 + t237) * t172) * MDP(16) - t271 * (qJD(2) * t158 + qJDD(2) * t170); qJDD(2) * MDP(2) + (t199 + t216) * MDP(3) + (-t212 * t250 + t214 * t228) * MDP(4) + (t155 * t159 - t156 * t163 + (t143 * t209 + t144 * t206 + t216) * pkin(2)) * MDP(5) + (t188 * t229 - t228 * t198 + t243 * t266 + t225 - t260) * MDP(8) + (t140 * t190 - t153 * t159 - g(3) * (pkin(2) * t214 + pkin(3) * t198 + qJ(4) * t196) + (t136 * t188 + t150 * t243) * t208 + (-t135 * t188 - t149 * t243) * t205 + t228 * (pkin(2) * t212 + pkin(3) * t196 - qJ(4) * t198)) * MDP(9) + (t145 * t173 - t164 * t165) * MDP(10) + (-t145 * t171 - t146 * t173 + t160 * t165 - t164 * t166) * MDP(11) + t147 * MDP(12) + t148 * MDP(13) + (t224 * qJDD(5) + t175 * t146 + t137 * t171 + t151 * t166 - t159 * t160 + t218 * t197 + (-qJD(5) * t223 - t243 * t173) * qJD(5)) * MDP(15) + (-t223 * qJDD(5) + t175 * t145 + t137 * t173 - t151 * t165 - t159 * t164 - t218 * t195 + (-qJD(5) * t224 + t243 * t171) * qJD(5)) * MDP(16) + t271 * (qJD(2) * t159 - qJDD(2) * t190 - t140 + t218); (qJDD(3) + t230) * MDP(5) + (t135 * t208 + t136 * t205 + t230) * MDP(9) + t148 * MDP(15) - t147 * MDP(16); (-t267 * qJD(2) - t218 + t221) * MDP(9) - t184 * MDP(15) + t235 * MDP(16) - t249 * MDP(8) * t215 + (-t251 - pkin(3) * MDP(9) + (t211 * MDP(15) + MDP(7)) * t205) * qJDD(2) + ((t164 - t219) * MDP(15) + (-t160 - t233) * MDP(16)) * qJD(5); t164 * t160 * MDP(10) + (-t160 ^ 2 + t164 ^ 2) * MDP(11) + t235 * MDP(12) + t226 * MDP(13) + qJDD(5) * MDP(14) + (-t211 * t134 + t213 * t133 - t151 * t164 - g(1) * (-t195 * t252 + t197 * t207) - g(2) * (-t195 * t253 - t197 * t210) + t195 * t260) * MDP(15) + (-t213 * t134 - t211 * t133 + t151 * t160 - g(1) * (-t195 * t207 - t197 * t252) - g(2) * (t195 * t210 - t197 * t253) + t197 * t260) * MDP(16) + ((t160 - t233) * MDP(12) + (t164 + t219) * MDP(13)) * qJD(5);];
tau = t1;
