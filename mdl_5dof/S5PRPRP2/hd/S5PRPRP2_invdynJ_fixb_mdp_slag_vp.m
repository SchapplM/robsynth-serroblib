% Calculate vector of inverse dynamics joint torques for
% S5PRPRP2
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRPRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:11
% EndTime: 2019-12-05 15:31:15
% DurationCPUTime: 1.46s
% Computational Cost: add. (725->213), mult. (1463->279), div. (0->0), fcn. (925->6), ass. (0->106)
t187 = cos(qJ(4));
t217 = qJD(2) * qJD(4);
t207 = t187 * t217;
t186 = sin(qJ(4));
t214 = qJDD(2) * t186;
t260 = t207 + t214;
t183 = sin(pkin(8));
t184 = cos(pkin(8));
t160 = -pkin(3) * t184 - pkin(6) * t183 - pkin(2);
t153 = qJD(2) * t160 + qJD(3);
t228 = qJD(2) * t184;
t158 = qJ(3) * t228 + qJD(1) * t183;
t229 = qJD(2) * t183;
t210 = qJ(5) * t229;
t136 = t158 * t187 + (t153 - t210) * t186;
t259 = qJD(4) * t136;
t180 = pkin(7) + qJ(2);
t176 = sin(t180);
t177 = cos(t180);
t257 = g(1) * t176 - g(2) * t177;
t218 = qJD(2) * qJD(3);
t219 = qJ(3) * qJDD(2);
t192 = t218 + t219;
t239 = t186 * t184;
t146 = t176 * t239 + t177 * t187;
t148 = t176 * t187 - t177 * t239;
t256 = -g(1) * t148 + g(2) * t146;
t255 = t260 * pkin(4) * t183 + qJDD(5);
t250 = qJ(5) * t183;
t249 = qJDD(2) * pkin(2);
t174 = t184 * qJDD(1);
t150 = t183 * t192 - t174;
t248 = t150 * t183;
t246 = (pkin(4) * t187 + pkin(3)) * t184;
t245 = t176 * t186;
t244 = t177 * t186;
t178 = t183 ^ 2;
t188 = qJD(2) ^ 2;
t243 = t178 * t188;
t242 = t183 * t187;
t241 = t184 * MDP(5);
t240 = t184 * t187;
t203 = t187 * t153 - t158 * t186;
t135 = -t187 * t210 + t203;
t167 = -qJD(4) + t228;
t132 = -pkin(4) * t167 + t135;
t238 = -t135 + t132;
t221 = t184 * qJD(3);
t224 = qJD(4) * t187;
t237 = t160 * t224 + t187 * t221;
t166 = qJ(3) * t240;
t236 = t186 * t160 + t166;
t235 = t177 * pkin(2) + t176 * qJ(3);
t234 = t184 ^ 2 + t178;
t181 = t186 ^ 2;
t182 = t187 ^ 2;
t233 = t181 - t182;
t232 = MDP(10) * t178;
t231 = MDP(12) * t183;
t230 = MDP(16) * t183;
t227 = qJD(2) * t186;
t226 = qJD(4) * t167;
t225 = qJD(4) * t186;
t223 = qJD(5) * t183;
t215 = qJDD(2) * t184;
t165 = -qJDD(4) + t215;
t222 = t165 * MDP(13);
t220 = -qJD(4) - t167;
t216 = qJD(2) * qJD(5);
t213 = qJDD(2) * t187;
t212 = qJ(3) * t239;
t211 = qJ(5) * t242;
t175 = t184 * qJD(1);
t157 = qJ(3) * t229 - t175;
t209 = t157 * t229;
t208 = t186 * t217;
t206 = -pkin(2) * t176 + t177 * qJ(3);
t205 = pkin(4) * t186 + qJ(3);
t202 = (-t181 - t182) * MDP(16);
t152 = qJDD(2) * t160 + qJDD(3);
t144 = t187 * t152;
t151 = qJDD(1) * t183 + t184 * t192;
t130 = -pkin(4) * t165 - t186 * t151 + t144 + (-qJ(5) * qJDD(2) - t216) * t242 - t259;
t201 = -t130 - t259;
t191 = -t187 * t151 - t186 * t152 - t153 * t224 + t158 * t225;
t131 = (-t260 * qJ(5) - t186 * t216) * t183 - t191;
t200 = qJD(4) * t132 - t131;
t199 = qJD(2) * t220;
t198 = t165 - t215;
t197 = t165 + t215;
t196 = g(1) * t177 + g(2) * t176;
t194 = t151 * t184 + t248;
t193 = t157 * t183 + t158 * t184;
t189 = -t167 * t228 + t226 - t243;
t173 = qJDD(3) - t249;
t159 = t183 * t184 * t208;
t156 = t187 * t160;
t149 = t177 * t240 + t245;
t147 = -t176 * t240 + t244;
t141 = t205 * t229 + qJD(5) - t175;
t139 = -t186 * t250 + t236;
t138 = -t211 + t156 + (-qJ(3) * t186 - pkin(4)) * t184;
t137 = t150 + t255;
t134 = -t186 * t221 - t187 * t223 + (-t166 + (-t160 + t250) * t186) * qJD(4);
t133 = -t186 * t223 + (-t211 - t212) * qJD(4) + t237;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (-t150 * t184 - g(3)) * MDP(8) + t159 * MDP(15) + (-t137 * t184 - g(3)) * MDP(17) + (t151 * MDP(8) + (MDP(14) * t198 - MDP(15) * t226 + MDP(17) * t201) * t186 + (t198 * MDP(15) + t131 * MDP(17) + ((t167 - t228) * MDP(14) - t132 * MDP(17)) * qJD(4)) * t187) * t183; qJDD(2) * MDP(2) + t257 * MDP(3) + t196 * MDP(4) + (t192 * t234 + t194 - t196) * MDP(7) + (-t173 * pkin(2) - g(1) * t206 - g(2) * t235 + qJ(3) * t194 + qJD(3) * t193) * MDP(8) + (qJDD(2) * t182 - 0.2e1 * t186 * t207) * t178 * MDP(9) + 0.2e1 * (-t186 * t213 + t217 * t233) * t232 + (t159 + (t167 * t225 - t187 * t197) * t183) * MDP(11) + (t197 * t186 + (t167 + t228) * t224) * t231 + t184 * t222 + (-g(1) * t147 - g(2) * t149 - t144 * t184 - t156 * t165 + ((qJD(2) * t178 + t167 * t184) * qJ(3) + t193) * t224 + (-(-qJD(4) * t160 - t221) * t167 - (-qJD(4) * t153 - t151) * t184 + t178 * t218 + t248 + (qJDD(2) * t178 + t165 * t184) * qJ(3)) * t186) * MDP(14) + ((-qJD(4) * t212 + t237) * t167 + t236 * t165 - t191 * t184 - g(1) * t146 - g(2) * t148 + (t150 * t187 - t157 * t225) * t183 + (t187 * t218 + (-t208 + t213) * qJ(3)) * t178) * MDP(15) + ((-qJDD(2) * t138 + (-qJD(4) * t139 - t134) * qJD(2) + t201) * t187 + (-qJDD(2) * t139 + (qJD(4) * t138 - t133) * qJD(2) + t200) * t186 + t257) * t230 + (t131 * t139 + t136 * t133 + t130 * t138 + t132 * t134 - g(1) * (pkin(4) * t244 - t176 * t246 + t206) - g(2) * (pkin(4) * t245 + t177 * t246 + t235) + (t137 * t205 + t141 * (pkin(4) * t224 + qJD(3)) - t257 * (-qJ(5) - pkin(6))) * t183) * MDP(17) + (-t183 * MDP(6) + t241) * (-t173 + t257 + t249); (-t158 * t228 + qJDD(3) - t209 - t257) * MDP(8) + (-t141 * t229 - t257) * MDP(17) - t234 * MDP(7) * t188 + (-t165 * MDP(14) + t189 * MDP(15) + (-t136 * t228 - t201) * MDP(17)) * t187 + (t189 * MDP(14) + t165 * MDP(15) + (t132 * t228 - t200) * MDP(17)) * t186 + (-t241 - pkin(2) * MDP(8) + (MDP(6) + t202) * t183) * qJDD(2); t187 * t186 * MDP(9) * t243 - t233 * t188 * t232 + (t186 * t199 + t213) * t183 * MDP(11) + (t187 * t199 - t214) * t231 - t222 + (t144 + (t158 * t220 - t209) * t187 + (g(3) * t183 + t153 * t220 - t151) * t186 + t256) * MDP(14) + (-t203 * t167 + g(1) * t149 - g(2) * t147 + (g(3) * t187 + t157 * t227) * t183 + t191) * MDP(15) + (-pkin(4) * t213 + (pkin(4) * qJD(4) - t238) * t227) * t230 + (t238 * t136 + (t130 + (-qJD(2) * t141 * t187 + g(3) * t186) * t183 + t256) * pkin(4)) * MDP(17); t202 * t243 + (g(3) * t184 - t174 + (t219 + (t132 * t187 + t136 * t186 + qJD(3)) * qJD(2) - t196) * t183 + t255) * MDP(17);];
tau = t1;
