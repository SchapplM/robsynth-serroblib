% Calculate vector of inverse dynamics joint torques for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:36:31
% EndTime: 2019-12-05 17:36:34
% DurationCPUTime: 1.53s
% Computational Cost: add. (891->224), mult. (1740->298), div. (0->0), fcn. (1071->10), ass. (0->113)
t190 = cos(qJ(4));
t223 = qJD(1) * qJD(4);
t212 = t190 * t223;
t188 = sin(qJ(4));
t218 = qJDD(1) * t188;
t262 = t212 + t218;
t186 = cos(pkin(7));
t172 = -pkin(1) * t186 - pkin(2);
t183 = sin(pkin(8));
t185 = cos(pkin(8));
t158 = -pkin(3) * t185 - pkin(6) * t183 + t172;
t148 = t158 * qJD(1) + qJD(3);
t184 = sin(pkin(7));
t171 = pkin(1) * t184 + qJ(3);
t163 = t171 * qJD(1);
t155 = qJD(2) * t183 + t163 * t185;
t235 = qJD(1) * t183;
t215 = qJ(5) * t235;
t135 = t155 * t190 + (t148 - t215) * t188;
t261 = qJD(4) * t135;
t224 = qJD(1) * qJD(3);
t159 = qJDD(1) * t171 + t224;
t180 = qJ(1) + pkin(7);
t176 = sin(t180);
t177 = cos(t180);
t245 = t185 * t188;
t149 = t176 * t245 + t177 * t190;
t151 = -t176 * t190 + t177 * t245;
t259 = -g(2) * t149 + g(3) * t151;
t257 = pkin(4) * t188;
t254 = qJ(5) * t183;
t174 = t185 * qJDD(2);
t251 = t159 * t183;
t146 = -t174 + t251;
t253 = t146 * t183;
t250 = (pkin(4) * t190 + pkin(3)) * t185;
t249 = t177 * t188;
t178 = t183 ^ 2;
t192 = qJD(1) ^ 2;
t248 = t178 * t192;
t247 = t183 * t190;
t246 = t185 * MDP(5);
t244 = t185 * t190;
t210 = t190 * t148 - t155 * t188;
t134 = -t190 * t215 + t210;
t234 = qJD(1) * t185;
t169 = -qJD(4) + t234;
t131 = -pkin(4) * t169 + t134;
t243 = -t134 + t131;
t228 = qJD(4) * t190;
t231 = qJD(3) * t190;
t242 = t158 * t228 + t185 * t231;
t160 = t171 * t244;
t241 = t188 * t158 + t160;
t240 = t185 ^ 2 + t178;
t181 = t188 ^ 2;
t182 = t190 ^ 2;
t239 = t181 - t182;
t238 = MDP(10) * t178;
t237 = MDP(12) * t183;
t236 = MDP(16) * t183;
t233 = qJD(1) * t188;
t232 = qJD(3) * t185;
t230 = qJD(4) * t169;
t229 = qJD(4) * t188;
t227 = qJD(5) * t183;
t219 = qJDD(1) * t185;
t168 = -qJDD(4) + t219;
t226 = t168 * MDP(13);
t225 = -qJD(4) - t169;
t222 = qJD(1) * qJD(5);
t220 = qJDD(1) * t172;
t217 = qJDD(1) * t190;
t216 = qJ(5) * t247;
t214 = t185 * t229;
t175 = t185 * qJD(2);
t154 = t163 * t183 - t175;
t213 = t154 * t235;
t189 = sin(qJ(1));
t211 = -pkin(1) * t189 + t177 * qJ(3);
t209 = (-t181 - t182) * MDP(16);
t145 = t158 * qJDD(1) + qJDD(3);
t142 = t190 * t145;
t147 = qJDD(2) * t183 + t159 * t185;
t129 = -pkin(4) * t168 - t188 * t147 + t142 + (-qJ(5) * qJDD(1) - t222) * t247 - t261;
t208 = -t129 - t261;
t196 = -t188 * t145 - t190 * t147 - t148 * t228 + t155 * t229;
t130 = (-t262 * qJ(5) - t188 * t222) * t183 - t196;
t207 = qJD(4) * t131 - t130;
t206 = qJD(1) * t225;
t205 = t168 - t219;
t204 = t168 + t219;
t203 = t262 * pkin(4) * t183 + qJDD(5) - t174;
t191 = cos(qJ(1));
t202 = -pkin(1) * t191 - pkin(2) * t177;
t201 = g(2) * t177 + g(3) * t176;
t200 = g(2) * t176 - g(3) * t177;
t199 = g(2) * t191 + g(3) * t189;
t198 = t147 * t185 + t253;
t197 = t154 * t183 + t155 * t185;
t194 = -t169 * t234 + t230 - t248;
t162 = qJDD(3) + t220;
t161 = t214 * t235;
t157 = t190 * t158;
t152 = -t176 * t188 - t177 * t244;
t150 = t176 * t244 - t249;
t139 = qJD(5) - t175 + (pkin(4) * t233 + t163) * t183;
t138 = -t188 * t254 + t241;
t137 = t203 + t251;
t136 = -t216 + t157 + (-t171 * t188 - pkin(4)) * t185;
t133 = -t188 * t232 - t190 * t227 + (-t160 + (-t158 + t254) * t188) * qJD(4);
t132 = -t188 * t227 + (-t171 * t245 - t216) * qJD(4) + t242;
t1 = [qJDD(1) * MDP(1) + t199 * MDP(2) + (-g(2) * t189 + g(3) * t191) * MDP(3) + (t199 + (t184 ^ 2 + t186 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t159 * t240 + t198 + t200) * MDP(7) + (t162 * t172 - g(2) * (-qJ(3) * t176 + t202) - g(3) * (-pkin(2) * t176 + t211) + t198 * t171 + t197 * qJD(3)) * MDP(8) + (qJDD(1) * t182 - 0.2e1 * t188 * t212) * t178 * MDP(9) + 0.2e1 * (-t188 * t217 + t239 * t223) * t238 + (t161 + (t169 * t229 - t204 * t190) * t183) * MDP(11) + (t204 * t188 + (t169 + t234) * t228) * t237 + t185 * t226 + (-g(2) * t152 + g(3) * t150 - t142 * t185 - t157 * t168 + ((qJD(1) * t178 + t169 * t185) * t171 + t197) * t228 + (-(-qJD(4) * t158 - t232) * t169 - (-qJD(4) * t148 - t147) * t185 + t178 * t224 + t253 + (qJDD(1) * t178 + t168 * t185) * t171) * t188) * MDP(14) + ((-t171 * t214 + t242) * t169 + t241 * t168 - t196 * t185 - g(2) * t151 - g(3) * t149 + (t146 * t190 - t154 * t229) * t183 + (t171 * t217 + (-t171 * t229 + t231) * qJD(1)) * t178) * MDP(15) + ((-qJDD(1) * t136 + (-qJD(4) * t138 - t133) * qJD(1) + t208) * t190 + (-qJDD(1) * t138 + (qJD(4) * t136 - t132) * qJD(1) + t207) * t188 + t201) * t236 + (t130 * t138 + t135 * t132 + t129 * t136 + t131 * t133 - g(2) * (-t177 * t250 + t202) - g(3) * (pkin(4) * t249 + t211) + (-g(2) * (-qJ(3) - t257) - g(3) * (-pkin(2) - t250)) * t176 + (t137 * (t171 + t257) + t139 * (pkin(4) * t228 + qJD(3)) - t201 * (-qJ(5) - pkin(6))) * t183) * MDP(17) + (-t183 * MDP(6) + t246) * (-t162 + t201 - t220); (qJDD(2) - g(1)) * MDP(4) + (-t146 * t185 - g(1)) * MDP(8) + t161 * MDP(15) + (-t137 * t185 - g(1)) * MDP(17) + (t147 * MDP(8) + (t205 * MDP(14) - MDP(15) * t230 + t208 * MDP(17)) * t188 + (t205 * MDP(15) + t130 * MDP(17) + ((t169 - t234) * MDP(14) - t131 * MDP(17)) * qJD(4)) * t190) * t183; (-t155 * t234 + qJDD(3) - t201 - t213) * MDP(8) + (-t139 * t235 - t201) * MDP(17) - t240 * MDP(7) * t192 + (-t168 * MDP(14) + t194 * MDP(15) + (-t135 * t234 - t208) * MDP(17)) * t190 + (t194 * MDP(14) + t168 * MDP(15) + (t131 * t234 - t207) * MDP(17)) * t188 + (-t246 + t172 * MDP(8) + (MDP(6) + t209) * t183) * qJDD(1); t190 * t188 * MDP(9) * t248 - t239 * t192 * t238 + (t188 * t206 + t217) * t183 * MDP(11) + (t190 * t206 - t218) * t237 - t226 + (t142 + (t225 * t155 - t213) * t190 + (g(1) * t183 + t225 * t148 - t147) * t188 + t259) * MDP(14) + (-t210 * t169 - g(2) * t150 - g(3) * t152 + (g(1) * t190 + t154 * t233) * t183 + t196) * MDP(15) + (-pkin(4) * t217 + (pkin(4) * qJD(4) - t243) * t233) * t236 + (t243 * t135 + (t129 + (-qJD(1) * t139 * t190 + g(1) * t188) * t183 + t259) * pkin(4)) * MDP(17); t209 * t248 + (g(1) * t185 + t203 + (t159 + (t131 * t190 + t135 * t188) * qJD(1) + t200) * t183) * MDP(17);];
tau = t1;
