% Calculate vector of inverse dynamics joint torques for
% S4RRPR4
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RRPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:36
% EndTime: 2019-12-31 17:02:38
% DurationCPUTime: 0.94s
% Computational Cost: add. (726->172), mult. (1070->221), div. (0->0), fcn. (669->12), ass. (0->94)
t191 = qJDD(1) + qJDD(2);
t195 = qJD(1) + qJD(2);
t200 = sin(qJ(2));
t234 = qJDD(1) * t200;
t203 = cos(qJ(2));
t235 = qJD(2) * t203;
t142 = t191 * qJ(3) + t195 * qJD(3) + (qJD(1) * t235 + t234) * pkin(1);
t197 = sin(pkin(7));
t198 = cos(pkin(7));
t238 = t197 ^ 2 + t198 ^ 2;
t223 = t238 * t142;
t252 = pkin(1) * qJD(1);
t232 = t200 * t252;
t253 = t203 * pkin(1);
t239 = -qJD(2) * t232 + qJDD(1) * t253;
t226 = qJDD(3) - t239;
t254 = t191 * pkin(2);
t147 = t226 - t254;
t196 = qJ(1) + qJ(2);
t188 = cos(t196);
t255 = g(2) * t188;
t257 = t147 + t255;
t199 = sin(qJ(4));
t202 = cos(qJ(4));
t156 = t202 * t197 + t199 * t198;
t150 = t156 * t195;
t187 = sin(t196);
t180 = g(1) * t187;
t256 = t180 + t239;
t152 = t156 * qJD(4);
t237 = qJD(1) * t203;
t218 = -pkin(1) * t237 + qJD(3);
t194 = pkin(7) + qJ(4);
t185 = sin(t194);
t250 = t187 * t185;
t186 = cos(t194);
t249 = t187 * t186;
t248 = t188 * t185;
t247 = t188 * t186;
t246 = t199 * t197;
t243 = t202 * t198;
t242 = t257 * t197;
t241 = t188 * pkin(2) + t187 * qJ(3);
t240 = g(1) * t188 + g(2) * t187;
t236 = qJD(2) * t200;
t233 = pkin(1) * t236;
t230 = t195 * t246;
t229 = t195 * t243;
t228 = qJD(4) * t229 + t156 * t191;
t177 = -t198 * pkin(3) - pkin(2);
t227 = t195 * t236;
t225 = -t187 * pkin(2) + t188 * qJ(3);
t224 = pkin(6) * t191 + t142;
t171 = pkin(1) * t235 + qJD(3);
t222 = t238 * t171;
t221 = t238 * t191;
t220 = t195 * t232;
t219 = -t240 + t223;
t167 = t191 * t243;
t217 = -t191 * t246 + t167;
t129 = -qJD(4) * t230 + t228;
t130 = t195 * t152 - t217;
t148 = -t229 + t230;
t155 = -t243 + t246;
t151 = t155 * qJD(4);
t216 = (-t129 * t155 - t156 * t130 + t151 * t148 - t150 * t152) * MDP(12) + (t129 * t156 - t150 * t151) * MDP(11) + (-t151 * qJD(4) + t156 * qJDD(4)) * MDP(13) + (-t152 * qJD(4) - t155 * qJDD(4)) * MDP(14) + t191 * MDP(4);
t176 = t200 * pkin(1) + qJ(3);
t153 = (-pkin(6) - t176) * t197;
t189 = t198 * pkin(6);
t154 = t198 * t176 + t189;
t215 = t202 * t153 - t199 * t154;
t214 = t199 * t153 + t202 * t154;
t163 = (-pkin(6) - qJ(3)) * t197;
t164 = t198 * qJ(3) + t189;
t213 = t202 * t163 - t199 * t164;
t212 = t199 * t163 + t202 * t164;
t211 = -t255 + t256;
t139 = t177 * t191 + t226;
t146 = t177 * t195 + t218;
t210 = -g(1) * t250 + g(2) * t248 + t139 * t156 - t146 * t151;
t209 = g(1) * t249 - g(2) * t247 + t139 * t155 + t146 * t152;
t208 = t220 - t255;
t182 = -pkin(2) - t253;
t207 = pkin(1) * t227 + t182 * t191;
t205 = t218 * t238;
t204 = cos(qJ(1));
t201 = sin(qJ(1));
t170 = t198 * t180;
t162 = t177 - t253;
t158 = t195 * qJ(3) + t232;
t157 = -t195 * pkin(2) + t218;
t134 = t224 * t198;
t133 = t224 * t197;
t1 = [qJDD(1) * MDP(1) + (g(1) * t201 - g(2) * t204) * MDP(2) + (g(1) * t204 + g(2) * t201) * MDP(3) + ((t191 * t203 - t227) * pkin(1) + t211) * MDP(5) + (((-qJDD(1) - t191) * t200 + (-qJD(1) - t195) * t235) * pkin(1) + t240) * MDP(6) + (t170 + (-t207 - t257) * t198) * MDP(7) + ((t207 - t180) * t197 + t242) * MDP(8) + (t176 * t221 + t195 * t222 + t219) * MDP(9) + (t147 * t182 + t157 * t233 - g(1) * (-t201 * pkin(1) + t225) - g(2) * (t204 * pkin(1) + t241) + t158 * t222 + t176 * t223) * MDP(10) + (t148 * t233 + t162 * t130 + t215 * qJDD(4) + (-t214 * qJD(4) - t156 * t171) * qJD(4) + t209) * MDP(16) + (t150 * t233 + t162 * t129 - t214 * qJDD(4) + (-t215 * qJD(4) + t155 * t171) * qJD(4) + t210) * MDP(17) + t216; (t208 + t256) * MDP(5) + ((-t234 + (-qJD(2) + t195) * t237) * pkin(1) + t240) * MDP(6) + (t170 + (-t147 + t208 + t254) * t198) * MDP(7) + ((-t220 - t254 - t180) * t197 + t242) * MDP(8) + (qJ(3) * t221 + t205 * t195 + t219) * MDP(9) + (-t147 * pkin(2) - g(1) * t225 - g(2) * t241 + qJ(3) * t223 - t157 * t232 + t205 * t158) * MDP(10) + (t177 * t130 + t213 * qJDD(4) + (-t156 * qJD(3) - t212 * qJD(4)) * qJD(4) + (-t200 * t148 + t203 * t152) * t252 + t209) * MDP(16) + (t177 * t129 - t212 * qJDD(4) + (t155 * qJD(3) - t213 * qJD(4)) * qJD(4) + (-t200 * t150 - t203 * t151) * t252 + t210) * MDP(17) + t216; (-t238 * t158 * t195 + qJDD(3) - t211) * MDP(10) - t167 * MDP(16) + t228 * MDP(17) + (-pkin(2) * MDP(10) - t198 * MDP(7) + (MDP(16) * t199 + MDP(8)) * t197) * t191 - t238 * MDP(9) * t195 ^ 2 + (0.2e1 * t150 * MDP(16) + (-t148 - t230) * MDP(17)) * qJD(4); t150 * t148 * MDP(11) + (-t148 ^ 2 + t150 ^ 2) * MDP(12) + t217 * MDP(14) + qJDD(4) * MDP(15) + (g(1) * t248 + g(2) * t250 - g(3) * t186 - t202 * t133 - t199 * t134 - t146 * t150) * MDP(16) + (g(1) * t247 + g(2) * t249 + g(3) * t185 + t199 * t133 - t202 * t134 + t146 * t148) * MDP(17) + (t228 + (t148 - t230) * qJD(4)) * MDP(13);];
tau = t1;
