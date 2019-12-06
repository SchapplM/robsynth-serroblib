% Calculate vector of inverse dynamics joint torques for
% S5PRPRP6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRPRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:25
% EndTime: 2019-12-05 15:41:28
% DurationCPUTime: 1.20s
% Computational Cost: add. (615->201), mult. (1081->244), div. (0->0), fcn. (602->6), ass. (0->107)
t171 = cos(qJ(4));
t169 = sin(qJ(4));
t194 = pkin(4) * t171 + qJ(5) * t169;
t137 = qJD(4) * t194 - qJD(5) * t171 + qJD(3);
t193 = pkin(4) * t169 - qJ(5) * t171;
t148 = qJ(3) + t193;
t210 = qJDD(2) * t148;
t170 = sin(qJ(2));
t213 = qJDD(1) * t170;
t172 = cos(qJ(2));
t217 = t172 * qJD(1);
t125 = t213 + t210 + (t137 + t217) * qJD(2);
t174 = qJD(4) ^ 2;
t167 = sin(pkin(7));
t168 = cos(pkin(7));
t195 = g(1) * t168 + g(2) * t167;
t216 = qJD(1) * qJD(2);
t246 = pkin(2) + pkin(6);
t177 = -t172 * (t195 + t216) - g(3) * t170 + t246 * t174;
t252 = qJD(2) * t137 + t125 + t177 + t210;
t196 = qJD(3) - t217;
t146 = -qJD(2) * t246 + t196;
t237 = t146 * t171;
t197 = qJD(5) - t237;
t243 = qJD(4) * pkin(4);
t132 = t197 - t243;
t220 = qJD(4) * qJ(5);
t238 = t146 * t169;
t134 = t220 + t238;
t160 = t170 * t216;
t212 = qJDD(1) * t172;
t189 = qJDD(3) + t160 - t212;
t250 = qJDD(2) * t246;
t135 = t189 - t250;
t130 = t169 * t135;
t207 = qJDD(4) * qJ(5);
t126 = t207 + t130 + (qJD(5) + t237) * qJD(4);
t131 = t171 * t135;
t219 = qJD(4) * t169;
t240 = qJDD(4) * pkin(4);
t249 = t146 * t219 - t240;
t127 = qJDD(5) - t131 + t249;
t191 = t126 * t169 - t127 * t171;
t178 = (t132 * t169 + t134 * t171) * qJD(4) + t191;
t164 = g(3) * t172;
t235 = t168 * t170;
t236 = t167 * t170;
t200 = -g(1) * t235 - g(2) * t236 + t164;
t251 = t178 + t200;
t247 = -(qJDD(1) - g(3)) * t170 + t172 * t195;
t242 = qJ(3) * t172;
t241 = qJDD(2) * pkin(2);
t239 = t134 * t169;
t234 = t169 * t170;
t233 = t169 * t171;
t232 = t170 * t171;
t165 = t169 ^ 2;
t166 = t171 ^ 2;
t230 = t165 - t166;
t229 = t165 + t166;
t175 = qJD(2) ^ 2;
t228 = t174 + t175;
t227 = qJD(1) * t170;
t226 = qJD(2) * qJ(3);
t223 = qJD(2) * t148;
t133 = t223 + t227;
t225 = qJD(2) * t133;
t153 = t226 + t227;
t222 = qJD(2) * t153;
t221 = qJD(2) * t171;
t218 = qJD(4) * t171;
t215 = qJD(2) * qJD(3);
t214 = qJD(2) * qJD(4);
t211 = qJDD(2) * qJ(3);
t209 = qJDD(2) * t170;
t208 = qJDD(2) * t171;
t206 = qJDD(4) * t169;
t205 = qJDD(4) * t246;
t204 = MDP(13) + MDP(15);
t203 = MDP(14) - MDP(17);
t202 = MDP(16) * qJDD(2);
t201 = t175 * t233;
t199 = g(3) * (t172 * pkin(2) + t170 * qJ(3));
t198 = t229 * MDP(16);
t141 = t167 * t171 + t168 * t234;
t143 = -t167 * t234 + t168 * t171;
t188 = g(1) * t141 - g(2) * t143 - t130;
t187 = -t200 + t212;
t186 = t153 + t226 - t227;
t140 = t167 * t169 - t168 * t232;
t142 = t167 * t232 + t168 * t169;
t185 = g(1) * t140 - g(2) * t142 + t171 * t164 + t131;
t184 = -qJDD(4) * t172 + 0.2e1 * t170 * t214;
t183 = t172 * t228 + t209;
t182 = qJDD(3) - t187;
t180 = -qJDD(5) + t185;
t179 = (t133 + t223 - t227) * qJD(4);
t136 = t211 + t213 + (qJD(3) + t217) * qJD(2);
t176 = t136 + t177 + t211 + t215;
t161 = qJDD(4) * t171;
t156 = t171 * t205;
t155 = t168 * t242;
t154 = t167 * t242;
t151 = -qJD(2) * pkin(2) + t196;
t147 = t194 * qJD(2);
t139 = t189 - t241;
t1 = [qJDD(1) * MDP(1) + (-MDP(3) + MDP(5)) * (-qJDD(2) * t172 + t170 * t175) + (-MDP(4) + MDP(6)) * (t172 * t175 + t209) + t203 * (-t169 * t184 + t171 * t183) + t204 * (t169 * t183 + t171 * t184) + (-MDP(1) - MDP(7) - MDP(18)) * g(3) + ((qJD(2) * t151 + t136) * MDP(7) + (qJD(2) * t239 - t132 * t221 + t125) * MDP(18) - t175 * t198) * t170 + ((-t139 + t222) * MDP(7) + (-t132 * t219 - t134 * t218 - t191 + t225) * MDP(18) + t229 * t202) * t172; qJDD(2) * MDP(2) + t187 * MDP(3) + t247 * MDP(4) + (t182 - 0.2e1 * t241) * MDP(5) + (0.2e1 * t211 + 0.2e1 * t215 - t247) * MDP(6) + (t136 * qJ(3) + t153 * qJD(3) - t139 * pkin(2) - g(1) * (-pkin(2) * t235 + t155) - g(2) * (-pkin(2) * t236 + t154) - t199 + (-t151 * t170 - t153 * t172) * qJD(1)) * MDP(7) + (qJDD(2) * t166 - 0.2e1 * t214 * t233) * MDP(8) + 0.2e1 * (-t169 * t208 + t214 * t230) * MDP(9) + (-t169 * t174 + t161) * MDP(10) + (-t171 * t174 - t206) * MDP(11) + (t169 * t176 + t186 * t218 - t156) * MDP(13) + ((-qJD(4) * t186 + t205) * t169 + t176 * t171) * MDP(14) + (t252 * t169 + t171 * t179 - t156) * MDP(15) + ((t160 + t250) * t229 - t251) * MDP(16) + ((t179 - t205) * t169 - t252 * t171) * MDP(17) + (t125 * t148 + t133 * t137 - g(1) * t155 - g(2) * t154 - t199 + (-g(3) * pkin(6) - t133 * qJD(1) - t195 * t193) * t172 - t178 * t246 + (-g(3) * t193 + (t132 * t171 - t239) * qJD(1) + t195 * t246) * t170) * MDP(18); -t175 * MDP(6) + (t160 + t182 - t222) * MDP(7) + (-t225 + t251) * MDP(18) + t204 * (-t169 * t228 + t161) - t203 * (t171 * t228 + t206) + (-pkin(2) * MDP(7) + MDP(5) - t198) * qJDD(2); MDP(8) * t201 - t230 * MDP(9) * t175 + MDP(10) * t208 - t169 * qJDD(2) * MDP(11) + qJDD(4) * MDP(12) + (-t153 * t221 + t185) * MDP(13) + ((t222 - t164) * t169 + t188) * MDP(14) + (0.2e1 * t240 + (-t133 * t171 - t147 * t169) * qJD(2) + t180) * MDP(15) + (-t194 * qJDD(2) + ((t134 - t220) * t171 + (-qJD(5) + t132 + t243) * t169) * qJD(2)) * MDP(16) + (t169 * t164 + 0.2e1 * t207 + 0.2e1 * qJD(4) * qJD(5) + (-t133 * t169 + t147 * t171) * qJD(2) - t188) * MDP(17) + (t126 * qJ(5) - t127 * pkin(4) - t133 * t147 - t132 * t238 - g(1) * (-pkin(4) * t140 + qJ(5) * t141) - g(2) * (pkin(4) * t142 - qJ(5) * t143) + t194 * t164 + t197 * t134) * MDP(18); (-qJDD(4) + t201) * MDP(15) + t171 * t202 + (-t166 * t175 - t174) * MDP(17) + (-qJD(4) * t134 + t133 * t221 - t180 + t249) * MDP(18);];
tau = t1;
