% Calculate vector of inverse dynamics joint torques for
% S4PRRR6
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4PRRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:09
% EndTime: 2019-12-31 16:35:11
% DurationCPUTime: 1.17s
% Computational Cost: add. (506->177), mult. (1103->253), div. (0->0), fcn. (785->10), ass. (0->92)
t180 = sin(qJ(2));
t176 = sin(pkin(7));
t177 = cos(pkin(7));
t203 = g(1) * t177 + g(2) * t176;
t197 = t203 * t180;
t183 = cos(qJ(2));
t235 = g(3) * t183;
t243 = t197 - t235;
t182 = cos(qJ(3));
t215 = qJD(2) * qJD(3);
t205 = t182 * t215;
t179 = sin(qJ(3));
t214 = qJDD(2) * t179;
t242 = t205 + t214;
t178 = sin(qJ(4));
t181 = cos(qJ(4));
t148 = t178 * t182 + t179 * t181;
t172 = qJD(3) + qJD(4);
t240 = t148 * t172;
t241 = qJD(2) * t240;
t213 = qJDD(2) * t182;
t126 = t178 * t214 - t181 * t213 + t241;
t147 = t178 * t179 - t181 * t182;
t194 = t147 * t172;
t220 = qJD(3) * t182;
t239 = -qJD(4) * t182 - t220;
t216 = qJD(1) * qJD(2);
t233 = qJDD(2) * pkin(2);
t144 = -qJDD(1) * t183 + t180 * t216 - t233;
t184 = qJD(3) ^ 2;
t238 = -pkin(5) * t184 + (t203 + t216) * t180 - t144 + t233 - t235;
t237 = pkin(5) + pkin(6);
t236 = g(3) * t180;
t234 = qJD(2) * pkin(2);
t157 = qJD(2) * pkin(5) + qJD(1) * t180;
t204 = pkin(6) * qJD(2) + t157;
t140 = t204 * t182;
t232 = t140 * t181;
t171 = qJDD(3) + qJDD(4);
t231 = t147 * t171;
t230 = t148 * t171;
t229 = t176 * t183;
t228 = t177 * t183;
t227 = qJDD(1) - g(3);
t173 = t179 ^ 2;
t226 = -t182 ^ 2 + t173;
t185 = qJD(2) ^ 2;
t225 = t184 + t185;
t224 = qJD(1) * t183;
t223 = qJD(2) * t179;
t222 = qJD(2) * t182;
t221 = qJD(3) * t179;
t219 = qJD(4) * t178;
t218 = qJD(4) * t181;
t212 = qJDD(2) * t183;
t211 = qJDD(3) * t179;
t210 = pkin(3) * t221;
t168 = -pkin(3) * t182 - pkin(2);
t209 = qJD(3) * t237;
t208 = t178 * t223;
t207 = t181 * t222;
t206 = t179 * t215;
t202 = g(1) * t176 - g(2) * t177;
t139 = t204 * t179;
t134 = qJD(3) * pkin(3) - t139;
t200 = -t134 * t178 - t232;
t151 = t237 * t179;
t152 = t237 * t182;
t199 = -t151 * t181 - t152 * t178;
t198 = -t151 * t178 + t152 * t181;
t196 = t203 * t183;
t193 = t206 - t213;
t125 = qJD(4) * t207 - t172 * t208 + t178 * t213 + t242 * t181;
t141 = -t207 + t208;
t143 = -t178 * t222 - t181 * t223;
t191 = -t143 * t141 * MDP(12) + (t141 * t172 + t125) * MDP(14) + (-t143 * t172 - t126) * MDP(15) + (-t141 ^ 2 + t143 ^ 2) * MDP(13) + t171 * MDP(16);
t158 = -t224 - t234;
t189 = -pkin(5) * qJDD(3) + (t158 + t224 - t234) * qJD(3);
t145 = qJDD(2) * pkin(5) + qJDD(1) * t180 + t183 * t216;
t188 = -t158 * qJD(2) - t145 + t196 + t236;
t128 = qJDD(3) * pkin(3) - t242 * pkin(6) - t179 * t145 - t157 * t220;
t146 = t168 * qJD(2) - t224;
t175 = qJ(3) + qJ(4);
t169 = sin(t175);
t170 = cos(t175);
t187 = -g(1) * (-t169 * t176 - t170 * t228) - g(2) * (t169 * t177 - t170 * t229) + t146 * t141 + t140 * t219 + t170 * t236 + (-t140 * t172 - t128) * t178;
t129 = -t193 * pkin(6) + t182 * t145 - t157 * t221;
t186 = -g(1) * (-t169 * t228 + t170 * t176) - g(2) * (-t169 * t229 - t170 * t177) + t200 * qJD(4) + t181 * t128 - t178 * t129 + t146 * t143 + t169 * t236;
t150 = t182 * t209;
t149 = t179 * t209;
t132 = t193 * pkin(3) + t144;
t1 = [t227 * MDP(1) + (-t180 * t185 + t212) * MDP(3) + (-qJDD(2) * t180 - t183 * t185) * MDP(4) + ((-0.2e1 * t206 + t213) * t183 + (-t225 * t182 - t211) * t180) * MDP(10) + ((-qJDD(3) * t180 - 0.2e1 * t183 * t215) * t182 + (t225 * t180 - t212) * t179) * MDP(11) + ((-t126 - t241) * t183 + ((t178 * t221 + t179 * t219 + t239 * t181) * t172 - t230 + qJD(2) * t141) * t180) * MDP(17) + ((qJD(2) * t194 - t125) * t183 + (-(t239 * t178 - t179 * t218 - t181 * t221) * t172 + t231 - qJD(2) * t143) * t180) * MDP(18); qJDD(2) * MDP(2) + (t227 * t183 + t197) * MDP(3) + (-t227 * t180 + t196) * MDP(4) + (qJDD(2) * t173 + 0.2e1 * t179 * t205) * MDP(5) + 0.2e1 * (t179 * t213 - t226 * t215) * MDP(6) + (t182 * t184 + t211) * MDP(7) + (qJDD(3) * t182 - t179 * t184) * MDP(8) + (t189 * t179 + t238 * t182) * MDP(10) + (-t238 * t179 + t189 * t182) * MDP(11) + (t125 * t148 + t143 * t194) * MDP(12) + (-t125 * t147 - t126 * t148 + t141 * t194 + t143 * t240) * MDP(13) + (-t172 * t194 + t230) * MDP(14) + (-t172 * t240 - t231) * MDP(15) + ((-t198 * qJD(4) + t149 * t178 - t150 * t181) * t172 + t199 * t171 + t141 * t210 + t168 * t126 + t132 * t147 + t146 * t240 + t243 * t170 + (-t180 * t141 + t183 * t240) * qJD(1)) * MDP(17) + (-(t199 * qJD(4) - t149 * t181 - t150 * t178) * t172 - t198 * t171 - t143 * t210 + t168 * t125 + t132 * t148 - t146 * t194 - t243 * t169 + (t180 * t143 - t183 * t194) * qJD(1)) * MDP(18); MDP(7) * t214 + MDP(8) * t213 + qJDD(3) * MDP(9) + (t188 * t179 - t202 * t182) * MDP(10) + (t202 * t179 + t188 * t182) * MDP(11) + (-(t139 * t178 - t232) * t172 + (-t141 * t223 + t181 * t171 - t172 * t219) * pkin(3) + t186) * MDP(17) + ((-qJD(4) * t134 - t139 * t172 - t129) * t181 + (t143 * t223 - t178 * t171 - t172 * t218) * pkin(3) + t187) * MDP(18) + t191 + (-t179 * t182 * MDP(5) + t226 * MDP(6)) * t185; (-t200 * t172 + t186) * MDP(17) + ((-t129 + (-qJD(4) + t172) * t134) * t181 + t187) * MDP(18) + t191;];
tau = t1;
