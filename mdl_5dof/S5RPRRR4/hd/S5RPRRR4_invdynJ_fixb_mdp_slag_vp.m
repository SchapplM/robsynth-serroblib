% Calculate vector of inverse dynamics joint torques for
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:41
% EndTime: 2019-12-05 18:14:43
% DurationCPUTime: 0.76s
% Computational Cost: add. (1092->154), mult. (1907->195), div. (0->0), fcn. (1095->14), ass. (0->96)
t188 = qJ(1) + pkin(9) + qJ(3);
t181 = qJ(4) + t188;
t175 = sin(t181);
t176 = cos(t181);
t260 = g(2) * t176 + g(3) * t175;
t194 = cos(pkin(9));
t180 = t194 * pkin(1) + pkin(2);
t165 = t180 * qJDD(1);
t197 = sin(qJ(3));
t201 = cos(qJ(3));
t193 = sin(pkin(9));
t254 = pkin(1) * t193;
t233 = qJD(1) * t254;
t225 = t197 * t233;
t167 = t180 * qJD(1);
t239 = qJD(3) * t167;
t142 = (qJDD(1) * t254 + t239) * t201 - qJD(3) * t225 + t197 * t165;
t149 = t201 * t167 - t225;
t190 = qJD(1) + qJD(3);
t146 = t190 * pkin(3) + t149;
t200 = cos(qJ(4));
t189 = qJDD(1) + qJDD(3);
t222 = t201 * t165 - t197 * t239;
t231 = qJD(1) * qJD(3) * t201;
t207 = (-qJDD(1) * t197 - t231) * t254 + t222;
t139 = t189 * pkin(3) + t207;
t150 = t197 * t167 + t201 * t233;
t196 = sin(qJ(4));
t238 = qJD(4) * t196;
t217 = -g(2) * t175 + g(3) * t176 - t196 * t139 + t150 * t238;
t209 = t217 - (qJD(4) * t146 + t142) * t200;
t168 = t201 * t180;
t257 = -t197 * t254 + t168;
t155 = pkin(3) + t257;
t156 = t197 * t180 + t201 * t254;
t241 = t196 * t155 + t200 * t156;
t178 = sin(t188);
t179 = cos(t188);
t258 = g(2) * t179 + g(3) * t178;
t245 = t196 * t150;
t141 = t200 * t149 - t245;
t250 = t196 * pkin(3);
t182 = pkin(8) + t250;
t249 = t200 * pkin(3);
t183 = -pkin(4) - t249;
t187 = qJD(4) + t190;
t256 = -qJDD(5) * t182 + (-qJD(4) * t249 + t183 * t187 + t141) * qJD(5);
t244 = t200 * t150;
t137 = t196 * t146 + t244;
t255 = t137 * qJD(4);
t253 = pkin(3) * t187;
t186 = qJDD(4) + t189;
t251 = t186 * pkin(4);
t153 = t257 * qJD(3);
t154 = t156 * qJD(3);
t248 = (t241 * qJD(4) + t196 * t153 + t200 * t154) * t187;
t247 = t137 * t187;
t246 = (t196 * t149 + t244) * t187;
t243 = qJDD(2) - g(1);
t229 = -t200 * t139 + t196 * t142;
t129 = t229 - t251 + t255;
t136 = t200 * t146 - t245;
t134 = -t187 * pkin(4) - t136;
t195 = sin(qJ(5));
t199 = cos(qJ(5));
t236 = t199 * qJD(5);
t242 = t129 * t195 + t134 * t236;
t191 = t195 ^ 2;
t240 = -t199 ^ 2 + t191;
t237 = qJD(5) * t187;
t232 = t134 * qJD(5) * t195 + t260 * t199;
t230 = -t146 - t253;
t198 = sin(qJ(1));
t202 = cos(qJ(1));
t223 = g(2) * t202 + g(3) * t198;
t203 = qJD(5) ^ 2;
t163 = qJDD(5) * t195 + t203 * t199;
t164 = qJDD(5) * t199 - t203 * t195;
t221 = 0.2e1 * (t195 * t186 * t199 - t240 * t237) * MDP(12) + (0.2e1 * t195 * t187 * t236 + t191 * t186) * MDP(11) + t163 * MDP(13) + t164 * MDP(14) + t186 * MDP(8);
t220 = t200 * t155 - t196 * t156;
t219 = -t229 + t260;
t216 = -pkin(8) * t203 + t247 + t251;
t143 = -pkin(4) - t220;
t144 = pkin(8) + t241;
t215 = -t143 * t186 - t144 * t203 - t248;
t214 = -t186 * pkin(8) - t134 * t187 + t209;
t130 = t220 * qJD(4) + t200 * t153 - t196 * t154;
t212 = -qJDD(5) * t144 + (t143 * t187 - t130) * qJD(5);
t211 = -pkin(4) * t237 - pkin(8) * qJDD(5) + qJD(5) * t136;
t210 = -t182 * t203 - t183 * t186 - t238 * t253 + t246;
t208 = t219 - t255;
t206 = t232 * MDP(16) + t242 * MDP(17) + t221;
t205 = -g(2) * t178 + g(3) * t179 - t142;
t185 = t187 ^ 2;
t184 = t189 * MDP(5);
t1 = [qJDD(1) * MDP(1) + t223 * MDP(2) + (-g(2) * t198 + g(3) * t202) * MDP(3) + (t223 + (t193 ^ 2 + t194 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + t184 + (-t154 * t190 + t168 * t189 + (-t231 + (-qJDD(1) - t189) * t197) * t254 + t222 + t258) * MDP(6) + (-t153 * t190 - t156 * t189 + t205) * MDP(7) + (t220 * t186 + t208 - t248) * MDP(9) + (-t130 * t187 - t241 * t186 + t209) * MDP(10) + (t212 * t195 + (-t129 + t215) * t199 + t232) * MDP(16) + (t212 * t199 + (-t215 - t260) * t195 + t242) * MDP(17) + t221; t164 * MDP(16) - t163 * MDP(17) + t243 * MDP(4); t184 + (t150 * t190 + t207 + t258) * MDP(6) + (t149 * t190 + t205) * MDP(7) + (t186 * t249 + t219 + t246) * MDP(9) + (t141 * t187 - t200 * t142 - t186 * t250 + t217) * MDP(10) + (t230 * MDP(9) * t196 + (t230 * MDP(10) - t150 * MDP(9)) * t200) * qJD(4) + ((-t129 + t210) * MDP(16) + t256 * MDP(17)) * t199 + ((-t210 - t260) * MDP(17) + t256 * MDP(16)) * t195 + t206; (t208 + t247) * MDP(9) + (t136 * t187 + t209) * MDP(10) + ((-t129 + t216) * MDP(16) + t211 * MDP(17)) * t199 + (t211 * MDP(16) + (-t216 - t260) * MDP(17)) * t195 + t206; qJDD(5) * MDP(15) + t240 * MDP(12) * t185 + (t186 * MDP(14) + t243 * MDP(16) + t214 * MDP(17)) * t199 + (-t185 * t199 * MDP(11) + t186 * MDP(13) + t214 * MDP(16) - t243 * MDP(17)) * t195;];
tau = t1;
