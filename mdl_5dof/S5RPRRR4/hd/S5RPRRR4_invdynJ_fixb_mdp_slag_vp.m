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
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:34:47
% EndTime: 2022-01-23 09:34:49
% DurationCPUTime: 0.79s
% Computational Cost: add. (1092->154), mult. (1907->195), div. (0->0), fcn. (1095->14), ass. (0->96)
t191 = qJDD(1) + qJDD(3);
t199 = sin(qJ(3));
t196 = cos(pkin(9));
t182 = pkin(1) * t196 + pkin(2);
t167 = t182 * qJDD(1);
t203 = cos(qJ(3));
t169 = t182 * qJD(1);
t241 = qJD(3) * t169;
t224 = t203 * t167 - t199 * t241;
t233 = qJD(1) * qJD(3) * t203;
t195 = sin(pkin(9));
t257 = pkin(1) * t195;
t209 = (-qJDD(1) * t199 - t233) * t257 + t224;
t141 = t191 * pkin(3) + t209;
t235 = qJD(1) * t257;
t226 = t199 * t235;
t144 = (qJDD(1) * t257 + t241) * t203 - qJD(3) * t226 + t199 * t167;
t190 = qJ(1) + pkin(9) + qJ(3);
t183 = qJ(4) + t190;
t178 = cos(t183);
t198 = sin(qJ(4));
t202 = cos(qJ(4));
t264 = -g(2) * t178 + t202 * t141 - t198 * t144;
t188 = qJDD(4) + t191;
t251 = t188 * pkin(4);
t151 = t203 * t169 - t226;
t192 = qJD(1) + qJD(3);
t148 = pkin(3) * t192 + t151;
t152 = t169 * t199 + t203 * t235;
t246 = t152 * t202;
t139 = t148 * t198 + t246;
t258 = t139 * qJD(4);
t263 = -t251 + t258 - t264;
t177 = sin(t183);
t240 = qJD(4) * t198;
t221 = g(1) * t178 + g(2) * t177 - t198 * t141 + t152 * t240;
t212 = t221 - (qJD(4) * t148 + t144) * t202;
t170 = t203 * t182;
t260 = -t199 * t257 + t170;
t157 = pkin(3) + t260;
t158 = t182 * t199 + t203 * t257;
t243 = t198 * t157 + t202 * t158;
t180 = sin(t190);
t181 = cos(t190);
t261 = g(1) * t180 - g(2) * t181;
t247 = t152 * t198;
t143 = t151 * t202 - t247;
t255 = pkin(3) * t198;
t184 = pkin(8) + t255;
t254 = pkin(3) * t202;
t185 = -pkin(4) - t254;
t189 = qJD(4) + t192;
t259 = -qJDD(5) * t184 + (-qJD(4) * t254 + t185 * t189 + t143) * qJD(5);
t256 = pkin(3) * t189;
t172 = g(1) * t177;
t155 = t260 * qJD(3);
t156 = t158 * qJD(3);
t250 = (t243 * qJD(4) + t198 * t155 + t202 * t156) * t189;
t249 = t139 * t189;
t248 = (t151 * t198 + t246) * t189;
t245 = qJDD(2) - g(3);
t138 = t148 * t202 - t247;
t136 = -pkin(4) * t189 - t138;
t197 = sin(qJ(5));
t201 = cos(qJ(5));
t244 = t136 * qJD(5) * t197 + t201 * t172;
t193 = t197 ^ 2;
t242 = -t201 ^ 2 + t193;
t239 = qJD(5) * t189;
t238 = qJD(5) * t201;
t234 = t136 * t238 + t263 * t197;
t232 = -t148 - t256;
t200 = sin(qJ(1));
t204 = cos(qJ(1));
t225 = g(1) * t200 - g(2) * t204;
t205 = qJD(5) ^ 2;
t165 = qJDD(5) * t197 + t201 * t205;
t166 = qJDD(5) * t201 - t197 * t205;
t223 = 0.2e1 * (t197 * t188 * t201 - t242 * t239) * MDP(12) + (0.2e1 * t189 * t197 * t238 + t188 * t193) * MDP(11) + t165 * MDP(13) + t166 * MDP(14) + t188 * MDP(8);
t222 = t157 * t202 - t158 * t198;
t219 = t172 + t264;
t218 = pkin(8) * t205 - t249 - t251;
t145 = -pkin(4) - t222;
t146 = pkin(8) + t243;
t217 = t145 * t188 + t146 * t205 + t250;
t216 = -t188 * pkin(8) - t136 * t189 + t212;
t132 = t222 * qJD(4) + t202 * t155 - t198 * t156;
t214 = -qJDD(5) * t146 + (t145 * t189 - t132) * qJD(5);
t213 = -pkin(4) * t239 - pkin(8) * qJDD(5) + qJD(5) * t138;
t211 = t184 * t205 + t185 * t188 + t240 * t256 - t248;
t210 = g(1) * t181 + g(2) * t180 - t144;
t208 = t219 - t258;
t207 = t244 * MDP(16) + t234 * MDP(17) + t223;
t187 = t189 ^ 2;
t186 = t191 * MDP(5);
t1 = [qJDD(1) * MDP(1) + t225 * MDP(2) + (g(1) * t204 + g(2) * t200) * MDP(3) + (t225 + (t195 ^ 2 + t196 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + t186 + (-t156 * t192 + t170 * t191 + (-t233 + (-qJDD(1) - t191) * t199) * t257 + t224 + t261) * MDP(6) + (-t155 * t192 - t158 * t191 + t210) * MDP(7) + (t222 * t188 + t208 - t250) * MDP(9) + (-t132 * t189 - t243 * t188 + t212) * MDP(10) + (t214 * t197 + (-t217 - t263) * t201 + t244) * MDP(16) + (t214 * t201 + (t217 - t172) * t197 + t234) * MDP(17) + t223; t166 * MDP(16) - t165 * MDP(17) + t245 * MDP(4); t186 + (t152 * t192 + t209 + t261) * MDP(6) + (t151 * t192 + t210) * MDP(7) + (t188 * t254 + t219 + t248) * MDP(9) + (t143 * t189 - t202 * t144 - t188 * t255 + t221) * MDP(10) + (t232 * MDP(9) * t198 + (t232 * MDP(10) - t152 * MDP(9)) * t202) * qJD(4) + ((t211 - t172) * MDP(17) + t259 * MDP(16)) * t197 + ((-t211 - t263) * MDP(16) + t259 * MDP(17)) * t201 + t207; (t208 + t249) * MDP(9) + (t138 * t189 + t212) * MDP(10) + (t213 * MDP(16) + (t218 - t172) * MDP(17)) * t197 + ((-t218 - t263) * MDP(16) + t213 * MDP(17)) * t201 + t207; qJDD(5) * MDP(15) + t242 * MDP(12) * t187 + (t188 * MDP(14) + t245 * MDP(16) + t216 * MDP(17)) * t201 + (-t187 * t201 * MDP(11) + t188 * MDP(13) + t216 * MDP(16) - t245 * MDP(17)) * t197;];
tau = t1;
