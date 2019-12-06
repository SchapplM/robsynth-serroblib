% Calculate vector of inverse dynamics joint torques for
% S5PRRPP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRPP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:59
% EndTime: 2019-12-05 16:07:03
% DurationCPUTime: 1.53s
% Computational Cost: add. (1169->238), mult. (2424->296), div. (0->0), fcn. (1547->8), ass. (0->100)
t263 = MDP(12) + MDP(15);
t202 = pkin(7) + qJ(2);
t196 = sin(t202);
t198 = cos(t202);
t261 = g(1) * t196 - g(2) * t198;
t228 = g(1) * t198 + g(2) * t196;
t207 = sin(pkin(8));
t208 = cos(pkin(8));
t211 = cos(qJ(3));
t245 = qJD(2) * t211;
t236 = t208 * t245;
t210 = sin(qJ(3));
t246 = qJD(2) * t210;
t172 = t207 * t246 - t236;
t179 = t207 * t211 + t208 * t210;
t175 = t179 * qJD(2);
t258 = pkin(3) * t211;
t195 = pkin(2) + t258;
t182 = -t195 * qJD(2) + qJD(4);
t147 = pkin(4) * t172 - qJ(5) * t175 + t182;
t203 = qJ(3) + pkin(8);
t197 = sin(t203);
t199 = cos(t203);
t200 = t211 * qJDD(1);
t253 = qJ(4) + pkin(6);
t230 = qJD(3) * t253;
t229 = qJD(2) * t230;
t259 = qJD(3) * qJD(1) + qJD(2) * qJD(4) + t253 * qJDD(2);
t146 = qJDD(3) * pkin(3) - t259 * t210 - t211 * t229 + t200;
t150 = (qJDD(1) - t229) * t210 + t259 * t211;
t134 = t208 * t146 - t207 * t150;
t232 = -qJDD(5) + t134;
t260 = -g(3) * t199 - t147 * t175 + t228 * t197 + t232;
t171 = t175 ^ 2;
t254 = g(3) * t211;
t252 = qJDD(3) * pkin(4);
t185 = t253 * t211;
t169 = qJD(1) * t210 + qJD(2) * t185;
t251 = t169 * t207;
t250 = t207 * t210;
t163 = t208 * t169;
t249 = t208 * t211;
t248 = qJDD(1) - g(3);
t135 = t207 * t146 + t208 * t150;
t235 = t253 * t210;
t167 = t211 * qJD(1) - qJD(2) * t235;
t166 = qJD(3) * pkin(3) + t167;
t149 = t207 * t166 + t163;
t205 = t210 ^ 2;
t247 = -t211 ^ 2 + t205;
t244 = qJD(3) * t210;
t154 = t167 * t208 - t251;
t243 = qJD(5) - t154;
t242 = qJD(2) * qJD(3);
t241 = qJDD(2) * t210;
t240 = qJDD(2) * t211;
t234 = t210 * t242;
t239 = pkin(3) * t234 + qJDD(4);
t238 = pkin(3) * t244;
t237 = qJDD(3) * qJ(5) + t135;
t233 = t211 * t242;
t226 = t195 * qJDD(2);
t225 = pkin(4) * t199 + qJ(5) * t197;
t148 = t166 * t208 - t251;
t223 = -0.2e1 * pkin(2) * t242 - pkin(6) * qJDD(3);
t174 = t179 * qJD(3);
t220 = -qJD(4) * t210 - t211 * t230;
t186 = t207 * t234;
t218 = t179 * qJDD(2) - t186;
t212 = qJD(3) ^ 2;
t217 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t212 + t261;
t213 = qJD(2) ^ 2;
t216 = pkin(2) * t213 - pkin(6) * qJDD(2) + t228;
t187 = t208 * t240;
t157 = qJD(2) * t174 + t207 * t241 - t187;
t158 = t208 * t233 + t218;
t215 = pkin(4) * t157 - qJ(5) * t158 - qJD(5) * t175 + t239;
t168 = qJD(4) * t211 - t210 * t230;
t153 = t168 * t207 - t208 * t220;
t155 = t208 * t168 + t207 * t220;
t160 = t185 * t207 + t208 * t235;
t161 = t208 * t185 - t207 * t235;
t214 = t153 * t175 - t155 * t172 - t161 * t157 + t158 * t160 - t228;
t194 = -pkin(3) * t208 - pkin(4);
t191 = pkin(3) * t207 + qJ(5);
t184 = qJDD(3) * t211 - t210 * t212;
t183 = qJDD(3) * t210 + t211 * t212;
t181 = t198 * t195;
t178 = -t249 + t250;
t177 = qJD(3) * t249 - t207 * t244;
t156 = pkin(4) * t178 - qJ(5) * t179 - t195;
t152 = t167 * t207 + t163;
t151 = pkin(3) * t246 + pkin(4) * t175 + qJ(5) * t172;
t143 = qJD(3) * qJ(5) + t149;
t141 = -qJD(3) * pkin(4) + qJD(5) - t148;
t138 = pkin(4) * t174 - qJ(5) * t177 - qJD(5) * t179 + t238;
t133 = -t232 - t252;
t132 = qJD(3) * qJD(5) + t237;
t131 = -t226 + t215;
t1 = [t248 * MDP(1) + t184 * MDP(10) - t183 * MDP(11) + (-t134 * t178 + t135 * t179 - t148 * t174 + t149 * t177 - g(3)) * MDP(13) + (-qJD(3) * t174 - qJDD(3) * t178) * MDP(14) + (qJD(3) * t177 + qJDD(3) * t179) * MDP(16) + (t132 * t179 + t133 * t178 + t141 * t174 + t143 * t177 - g(3)) * MDP(17) + t263 * (-t179 * t157 + t158 * t178 - t177 * t172 + t174 * t175); qJDD(2) * MDP(2) + t261 * MDP(3) + t228 * MDP(4) + (qJDD(2) * t205 + 0.2e1 * t210 * t233) * MDP(5) + 0.2e1 * (t210 * t240 - t247 * t242) * MDP(6) + t183 * MDP(7) + t184 * MDP(8) + (t223 * t210 + t217 * t211) * MDP(10) + (-t217 * t210 + t223 * t211) * MDP(11) + (-t134 * t179 - t135 * t178 - t148 * t177 - t149 * t174 + t214) * MDP(12) + (t135 * t161 + t149 * t155 - t134 * t160 - t148 * t153 - (-t226 + t239) * t195 + t182 * t238 - g(1) * (-t195 * t196 + t198 * t253) - g(2) * (t196 * t253 + t181)) * MDP(13) + (-qJD(3) * t153 - qJDD(3) * t160 + t131 * t178 + t138 * t172 + t147 * t174 + t156 * t157 + t199 * t261) * MDP(14) + (-t132 * t178 + t133 * t179 + t141 * t177 - t143 * t174 + t214) * MDP(15) + (qJD(3) * t155 + qJDD(3) * t161 - t131 * t179 - t138 * t175 - t147 * t177 - t156 * t158 + t197 * t261) * MDP(16) + (-g(2) * t181 + t131 * t156 + t132 * t161 + t133 * t160 + t147 * t138 + t141 * t153 + t143 * t155 + (-g(1) * t253 - g(2) * t225) * t198 + (-g(1) * (-t195 - t225) - g(2) * t253) * t196) * MDP(17); MDP(7) * t241 + MDP(8) * t240 + qJDD(3) * MDP(9) + (t216 * t210 + t200 - t254) * MDP(10) + (-t248 * t210 + t216 * t211) * MDP(11) + ((t149 - t152) * t175 + (-t148 + t154) * t172 + (-t157 * t207 - t158 * t208) * pkin(3)) * MDP(12) + (t148 * t152 - t149 * t154 + (-t254 + t134 * t208 + t135 * t207 + (-qJD(2) * t182 + t228) * t210) * pkin(3)) * MDP(13) + (qJD(3) * t152 - t151 * t172 + (pkin(4) - t194) * qJDD(3) + t260) * MDP(14) + (-t157 * t191 + t158 * t194 + (t143 - t152) * t175 + (t141 - t243) * t172) * MDP(15) + (-g(3) * t197 + qJDD(3) * t191 - t147 * t172 + t151 * t175 - t228 * t199 + (0.2e1 * qJD(5) - t154) * qJD(3) + t237) * MDP(16) + (t132 * t191 + t133 * t194 - t147 * t151 - t141 * t152 - g(3) * (t225 + t258) + t243 * t143 + t228 * (pkin(3) * t210 + pkin(4) * t197 - qJ(5) * t199)) * MDP(17) + (-t210 * t211 * MDP(5) + t247 * MDP(6)) * t213; (t148 * t175 + t149 * t172 + t239 - t261) * MDP(13) - t187 * MDP(14) + t186 * MDP(16) + (-t141 * t175 + t143 * t172 + t215 - t261) * MDP(17) + (MDP(14) * t250 - t179 * MDP(16) + (-MDP(13) - MDP(17)) * t195) * qJDD(2) + ((t207 * t245 + t208 * t246 + t175) * MDP(14) + (t172 - t236) * MDP(16)) * qJD(3) + t263 * (-t172 ^ 2 - t171); (t172 * t175 - qJDD(3)) * MDP(14) + ((t172 + t236) * qJD(3) + t218) * MDP(15) + (-t171 - t212) * MDP(16) + (-qJD(3) * t143 - t252 - t260) * MDP(17);];
tau = t1;
