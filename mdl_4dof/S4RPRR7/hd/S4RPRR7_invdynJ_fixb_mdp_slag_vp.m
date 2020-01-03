% Calculate vector of inverse dynamics joint torques for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RPRR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:08
% EndTime: 2019-12-31 16:54:11
% DurationCPUTime: 1.70s
% Computational Cost: add. (1014->226), mult. (2439->305), div. (0->0), fcn. (1824->10), ass. (0->107)
t196 = sin(pkin(7));
t199 = sin(qJ(3));
t236 = qJD(1) * t199;
t226 = t196 * t236;
t197 = cos(pkin(7));
t202 = cos(qJ(3));
t235 = qJD(1) * t202;
t183 = t197 * t235;
t229 = qJDD(1) * t202;
t230 = qJDD(1) * t199;
t227 = qJD(3) * t183 + t196 * t229 + t197 * t230;
t151 = -qJD(3) * t226 + t227;
t260 = qJD(3) * qJD(4) + t151;
t251 = pkin(5) + qJ(2);
t177 = t251 * t196;
t173 = qJD(1) * t177;
t178 = t251 * t197;
t174 = qJD(1) * t178;
t154 = -t173 * t199 + t174 * t202;
t232 = qJD(1) * qJD(2);
t255 = t251 * qJDD(1) + t232;
t160 = t255 * t196;
t161 = t255 * t197;
t216 = t160 * t202 + t161 * t199;
t135 = -qJDD(3) * pkin(3) + qJD(3) * t154 + t216;
t167 = t183 - t226;
t162 = qJD(4) - t167;
t172 = t196 * t202 + t197 * t199;
t168 = t172 * qJD(1);
t195 = pkin(7) + qJ(3);
t188 = sin(t195);
t189 = cos(t195);
t200 = sin(qJ(1));
t203 = cos(qJ(1));
t219 = g(1) * t203 + g(2) * t200;
t208 = -g(3) * t189 + t219 * t188;
t259 = t208 - (pkin(3) * t168 + t162 * pkin(6)) * t162 - t135;
t258 = qJ(2) * qJDD(1);
t218 = g(1) * t200 - g(2) * t203;
t257 = qJDD(2) - t218;
t171 = t196 * t199 - t202 * t197;
t214 = -t177 * t202 - t178 * t199;
t140 = -t171 * qJD(2) + t214 * qJD(3);
t170 = t172 * qJD(3);
t182 = t197 * t229;
t217 = -t196 * t230 + t182;
t152 = qJD(1) * t170 - t217;
t145 = qJDD(4) + t152;
t153 = -t173 * t202 - t174 * t199;
t146 = -qJD(3) * pkin(3) - t153;
t185 = -pkin(2) * t197 - pkin(1);
t150 = pkin(3) * t171 - pkin(6) * t172 + t185;
t156 = -t177 * t199 + t178 * t202;
t169 = t171 * qJD(3);
t215 = -t160 * t199 + t161 * t202;
t134 = qJDD(3) * pkin(6) + qJD(3) * t153 + t215;
t176 = t185 * qJD(1) + qJD(2);
t139 = -pkin(3) * t167 - pkin(6) * t168 + t176;
t223 = qJD(4) * t139 + t134;
t254 = t135 * t172 - t156 * t145 - t146 * t169 - (qJD(4) * t150 + t140) * t162 - t223 * t171;
t253 = g(3) * t188;
t250 = qJDD(1) * pkin(1);
t198 = sin(qJ(4));
t201 = cos(qJ(4));
t228 = t198 * qJDD(3) + t260 * t201;
t233 = qJD(4) * t198;
t137 = -t168 * t233 + t228;
t249 = t137 * t198;
t248 = t145 * t198;
t247 = t150 * t145;
t157 = -t201 * qJD(3) + t168 * t198;
t246 = t157 * t162;
t245 = t157 * t168;
t159 = qJD(3) * t198 + t168 * t201;
t244 = t159 * t162;
t243 = t159 * t168;
t242 = t197 * MDP(4);
t241 = t198 * t200;
t240 = t198 * t203;
t239 = t200 * t201;
t142 = t201 * t145;
t238 = t201 * t203;
t237 = t196 ^ 2 + t197 ^ 2;
t234 = qJD(4) * t172;
t224 = t162 * t201;
t175 = t185 * qJDD(1) + qJDD(2);
t136 = pkin(3) * t152 - pkin(6) * t151 + t175;
t147 = qJD(3) * pkin(6) + t154;
t222 = qJD(4) * t147 - t136;
t220 = 0.2e1 * t237;
t213 = t250 - t257;
t212 = t142 + (t167 * t198 - t233) * t162;
t211 = -t169 * t201 - t172 * t233;
t207 = -pkin(6) * t145 + (t146 + t153) * t162;
t206 = t220 * t232 - t219;
t191 = t201 * qJDD(3);
t166 = t189 * t238 + t241;
t165 = -t189 * t240 + t239;
t164 = -t189 * t239 + t240;
t163 = t189 * t241 + t238;
t149 = pkin(3) * t170 + pkin(6) * t169;
t141 = t172 * qJD(2) + t156 * qJD(3);
t138 = t159 * qJD(4) + t151 * t198 - t191;
t133 = t201 * t136;
t132 = t139 * t198 + t147 * t201;
t131 = t139 * t201 - t147 * t198;
t1 = [qJDD(1) * MDP(1) + t218 * MDP(2) + t219 * MDP(3) + (t220 * t258 + t206) * MDP(6) + (t213 * pkin(1) + (t237 * t258 + t206) * qJ(2)) * MDP(7) + (t151 * t172 - t168 * t169) * MDP(8) + (-t151 * t171 - t152 * t172 - t167 * t169 - t168 * t170) * MDP(9) + (-qJD(3) * t169 + qJDD(3) * t172) * MDP(10) + (-qJD(3) * t170 - qJDD(3) * t171) * MDP(11) + (-qJD(3) * t141 + qJDD(3) * t214 + t152 * t185 + t170 * t176 + t171 * t175 + t218 * t189) * MDP(13) + (-qJD(3) * t140 - qJDD(3) * t156 + t151 * t185 - t169 * t176 + t172 * t175 - t218 * t188) * MDP(14) + (t137 * t172 * t201 + t211 * t159) * MDP(15) + (-(-t157 * t201 - t159 * t198) * t169 + (-t249 - t138 * t201 + (t157 * t198 - t159 * t201) * qJD(4)) * t172) * MDP(16) + (t137 * t171 + t172 * t142 + t159 * t170 + t211 * t162) * MDP(17) + (-t172 * t248 - t138 * t171 - t157 * t170 + (t169 * t198 - t201 * t234) * t162) * MDP(18) + (t145 * t171 + t162 * t170) * MDP(19) + (-g(1) * t164 - g(2) * t166 + t131 * t170 + t133 * t171 - t214 * t138 + t141 * t157 + (t247 + t149 * t162 + (t146 * t172 - t147 * t171 - t156 * t162) * qJD(4)) * t201 + t254 * t198) * MDP(20) + (-g(1) * t163 - g(2) * t165 - t132 * t170 - t214 * t137 + t141 * t159 + (-(-qJD(4) * t156 + t149) * t162 - t247 + t222 * t171 - t146 * t234) * t198 + t254 * t201) * MDP(21) + (-MDP(5) * t196 + t242) * (t213 + t250); t257 * MDP(7) - t182 * MDP(13) + t227 * MDP(14) + (t212 - t245) * MDP(20) + (-t162 ^ 2 * t201 - t243 - t248) * MDP(21) + (-t242 - pkin(1) * MDP(7) + (MDP(13) * t199 + MDP(5)) * t196) * qJDD(1) + ((t196 * t235 + t197 * t236 + t168) * MDP(13) + (t167 - t226) * MDP(14)) * qJD(3) + (-MDP(7) * qJ(2) - MDP(6)) * qJD(1) ^ 2 * t237; -t167 ^ 2 * MDP(9) + ((-t167 - t226) * qJD(3) + t227) * MDP(10) + t217 * MDP(11) + qJDD(3) * MDP(12) + (t208 - t216) * MDP(13) + (-t167 * t176 + t189 * t219 - t215 + t253) * MDP(14) + (t159 * t224 + t249) * MDP(15) + ((t137 - t246) * t201 + (-t138 - t244) * t198) * MDP(16) + (t162 * t224 - t243 + t248) * MDP(17) + (t212 + t245) * MDP(18) + (-pkin(3) * t138 - t154 * t157 + t207 * t198 + t201 * t259) * MDP(20) + (-pkin(3) * t137 - t154 * t159 - t198 * t259 + t207 * t201) * MDP(21) + (-MDP(13) * t176 - MDP(19) * t162 - MDP(20) * t131 + MDP(21) * t132 - MDP(8) * t167 + MDP(9) * t168) * t168; t159 * t157 * MDP(15) + (-t157 ^ 2 + t159 ^ 2) * MDP(16) + (t228 + t246) * MDP(17) + (t191 + t244) * MDP(18) + t145 * MDP(19) + (-g(1) * t165 + g(2) * t163 + t132 * t162 - t146 * t159 + t133) * MDP(20) + (g(1) * t166 - g(2) * t164 + t131 * t162 + t146 * t157) * MDP(21) + ((-t134 + t253) * MDP(21) + (-MDP(18) * t168 - MDP(20) * t147 - MDP(21) * t139) * qJD(4)) * t201 + (-qJD(4) * t168 * MDP(17) - t260 * MDP(18) + (-t223 + t253) * MDP(20) + t222 * MDP(21)) * t198;];
tau = t1;
