% Calculate vector of inverse dynamics joint torques for
% S5RPPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPPPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:34
% EndTime: 2019-12-31 17:46:36
% DurationCPUTime: 1.36s
% Computational Cost: add. (661->198), mult. (1154->253), div. (0->0), fcn. (766->10), ass. (0->99)
t206 = sin(pkin(8));
t208 = cos(pkin(8));
t255 = t206 ^ 2 + t208 ^ 2;
t269 = t208 * MDP(10) - t206 * MDP(11);
t207 = sin(pkin(7));
t209 = cos(pkin(7));
t265 = sin(qJ(1));
t266 = cos(qJ(1));
t162 = -t265 * t207 - t266 * t209;
t164 = t266 * t207 - t265 * t209;
t229 = -g(1) * t162 - g(2) * t164;
t268 = qJDD(1) * pkin(3) + qJDD(4);
t267 = qJD(5) ^ 2;
t213 = -pkin(1) - pkin(2);
t264 = pkin(4) * t208;
t169 = t209 * qJ(2) + t207 * t213;
t165 = -qJ(4) + t169;
t261 = pkin(6) - t165;
t260 = pkin(1) * qJDD(1);
t214 = qJD(1) ^ 2;
t259 = t209 * t214;
t247 = qJD(1) * qJD(2);
t180 = t209 * t247;
t172 = t213 * qJDD(1) + qJDD(2);
t243 = qJDD(1) * t209;
t258 = qJ(2) * t243 + t207 * t172;
t148 = t180 + t258;
t142 = -qJDD(1) * qJ(4) - qJD(1) * qJD(4) + t148;
t135 = t206 * qJDD(3) + t208 * t142;
t173 = t213 * qJD(1) + qJD(2);
t254 = qJ(2) * qJD(1);
t155 = t207 * t173 + t209 * t254;
t257 = t266 * pkin(1) + t265 * qJ(2);
t256 = g(1) * t265 - g(2) * t266;
t211 = sin(qJ(5));
t253 = qJD(1) * t211;
t212 = cos(qJ(5));
t252 = qJD(1) * t212;
t251 = qJD(2) * t207;
t248 = qJ(2) * qJDD(1);
t246 = qJD(1) * qJD(5);
t245 = qJDD(1) * t207;
t244 = qJDD(1) * t208;
t242 = qJDD(1) * t211;
t241 = qJDD(1) * t212;
t161 = t206 * t211 - t212 * t208;
t240 = qJDD(5) * t161;
t163 = t206 * t212 + t208 * t211;
t239 = qJDD(5) * t163;
t238 = 0.2e1 * t247;
t237 = t266 * pkin(2) + t257;
t236 = t206 * t253;
t235 = t208 * t252;
t178 = t207 * t247;
t234 = -qJ(2) * t245 + t172 * t209;
t154 = t173 * t209 - t207 * t254;
t168 = -t207 * qJ(2) + t209 * t213;
t233 = qJDD(1) * t255;
t232 = qJDD(2) - t260;
t166 = pkin(3) - t168;
t231 = -t265 * pkin(1) + t266 * qJ(2);
t230 = g(1) * t164 - g(2) * t162;
t147 = -t178 + t234;
t174 = t208 * t241;
t228 = t206 * t242 - t174;
t189 = t208 * qJDD(3);
t134 = -t142 * t206 + t189;
t227 = -t134 * t206 + t135 * t208;
t226 = t255 * (-qJD(1) * qJ(4) + t155);
t152 = t261 * t206;
t153 = t261 * t208;
t225 = t152 * t212 + t153 * t211;
t224 = t152 * t211 - t153 * t212;
t223 = t154 * t207 - t155 * t209;
t159 = t161 * qJD(5);
t221 = qJD(5) * t159 - t239;
t160 = t163 * qJD(5);
t220 = qJD(5) * t160 + t240;
t150 = qJD(1) * pkin(3) + qJD(4) - t154;
t219 = g(1) * t266 + g(2) * t265;
t145 = -t147 + t268;
t218 = -t206 * t252 - t208 * t253;
t217 = -t265 * pkin(2) + t231;
t216 = -t230 - t234;
t204 = pkin(8) + qJ(5);
t193 = cos(t204);
t192 = sin(t204);
t176 = qJD(2) * t209 - qJD(4);
t171 = qJD(5) * t236;
t158 = t163 * qJD(1);
t157 = -t235 + t236;
t156 = t166 + t264;
t146 = qJD(1) * t264 + t150;
t141 = qJD(1) * t160 + t228;
t140 = -qJD(5) * t235 - t163 * qJDD(1) + t171;
t138 = pkin(4) * t244 + t145;
t133 = -pkin(6) * t244 + t135;
t132 = t189 + (pkin(6) * qJDD(1) - t142) * t206;
t1 = [qJDD(1) * MDP(1) + t256 * MDP(2) + t219 * MDP(3) + (-qJDD(2) + t256 + 0.2e1 * t260) * MDP(4) + (-t219 + t238 + 0.2e1 * t248) * MDP(5) + (-t232 * pkin(1) - g(1) * t231 - g(2) * t257 + (t238 + t248) * qJ(2)) * MDP(6) + (-qJDD(1) * t168 + 0.2e1 * t178 + t216) * MDP(7) + (qJDD(1) * t169 + 0.2e1 * t180 - t229 + t258) * MDP(8) + (-g(1) * t217 - g(2) * t237 - t223 * qJD(2) + t147 * t168 + t148 * t169) * MDP(9) + (-t255 * t176 * qJD(1) - t165 * t233 - t227 + t229) * MDP(12) + (t145 * t166 + t150 * t251 - g(1) * (t164 * pkin(3) + t162 * qJ(4) + t217) - g(2) * (-pkin(3) * t162 + qJ(4) * t164 + t237) + t226 * t176 + t227 * t165) * MDP(13) + (-t140 * t163 - t158 * t159) * MDP(14) + (t140 * t161 - t141 * t163 + t157 * t159 - t158 * t160) * MDP(15) + t221 * MDP(16) + t220 * MDP(17) + (-t157 * t251 - t156 * t141 - t138 * t161 - t146 * t160 + t225 * qJDD(5) - t230 * t193 + (-t224 * qJD(5) - t163 * t176) * qJD(5)) * MDP(19) + (-t158 * t251 + t156 * t140 - t138 * t163 + t146 * t159 - t224 * qJDD(5) + t230 * t192 + (-t225 * qJD(5) + t161 * t176) * qJD(5)) * MDP(20) + t269 * (qJDD(1) * t166 + t145 + t178 - t230); -qJDD(1) * MDP(4) - t214 * MDP(5) + (-qJ(2) * t214 + t232 - t256) * MDP(6) + (t245 - t259) * MDP(8) + (t223 * qJD(1) + t147 * t209 + t148 * t207 - t256) * MDP(9) + (-t207 * t233 + t255 * t259) * MDP(12) + (-t145 * t209 + t227 * t207 + (-t150 * t207 - t226 * t209) * qJD(1) - t256) * MDP(13) + ((t163 * t246 + t141) * t209 + (qJD(1) * t157 + t161 * t267 - t239) * t207) * MDP(19) + ((-t161 * t246 - t140) * t209 + (qJD(1) * t158 + t163 * t267 + t240) * t207) * MDP(20) + (-MDP(7) - t269) * (t207 * t214 + t243); (qJDD(3) + g(3)) * MDP(9) + (t134 * t208 + t135 * t206 + g(3)) * MDP(13) - t220 * MDP(19) + t221 * MDP(20); (t226 * qJD(1) + t178 + t216 + t268) * MDP(13) + t174 * MDP(19) + t171 * MDP(20) - t255 * MDP(12) * t214 + ((-t211 * MDP(20) + MDP(10)) * t208 + (-t211 * MDP(19) - t212 * MDP(20) - MDP(11)) * t206) * qJDD(1) + ((-t158 + t218) * MDP(19) + (t157 - t235) * MDP(20)) * qJD(5); t158 * t157 * MDP(14) + (-t157 ^ 2 + t158 ^ 2) * MDP(15) + (-t206 * t241 - t208 * t242 + t171) * MDP(16) + t228 * MDP(17) + qJDD(5) * MDP(18) + (g(3) * t193 + t212 * t132 - t211 * t133 + t146 * t158 + t192 * t229) * MDP(19) + (-g(3) * t192 - t211 * t132 - t212 * t133 - t146 * t157 + t193 * t229) * MDP(20) + ((-t157 - t235) * MDP(16) + (-t158 - t218) * MDP(17)) * qJD(5);];
tau = t1;
