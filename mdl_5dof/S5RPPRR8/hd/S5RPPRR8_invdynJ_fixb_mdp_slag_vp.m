% Calculate vector of inverse dynamics joint torques for
% S5RPPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPRR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:17
% EndTime: 2019-12-31 18:01:19
% DurationCPUTime: 1.01s
% Computational Cost: add. (830->157), mult. (1213->204), div. (0->0), fcn. (714->10), ass. (0->83)
t246 = qJD(1) - qJD(4);
t180 = t246 ^ 2;
t187 = sin(qJ(5));
t247 = 0.2e1 * t187;
t190 = cos(qJ(5));
t224 = qJD(5) * t190;
t244 = MDP(18) * t190 - MDP(19) * t187 + MDP(11);
t186 = cos(pkin(8));
t193 = -pkin(1) - pkin(2);
t185 = sin(pkin(8));
t237 = qJ(2) * t185;
t163 = t186 * t193 - t237;
t160 = -pkin(3) + t163;
t164 = qJ(2) * t186 + t185 * t193;
t188 = sin(qJ(4));
t191 = cos(qJ(4));
t232 = t188 * t160 + t191 * t164;
t168 = t193 * qJDD(1) + qJDD(2);
t162 = t186 * t168;
t216 = -pkin(3) - t237;
t222 = qJD(1) * qJD(2);
t218 = t185 * t222;
t137 = t216 * qJDD(1) + t162 - t218;
t217 = t186 * t222;
t221 = qJDD(1) * t186;
t231 = qJ(2) * t221 + t185 * t168;
t142 = t217 + t231;
t172 = t193 * qJD(1) + qJD(2);
t165 = t186 * t172;
t144 = t216 * qJD(1) + t165;
t227 = qJ(2) * qJD(1);
t147 = t172 * t185 + t186 * t227;
t189 = sin(qJ(1));
t192 = cos(qJ(1));
t220 = pkin(8) + qJ(4);
t214 = sin(t220);
t215 = cos(t220);
t148 = -t189 * t214 - t192 * t215;
t149 = -t189 * t215 + t192 * t214;
t225 = qJD(4) * t191;
t226 = qJD(4) * t188;
t199 = g(1) * t148 + g(2) * t149 + t188 * t137 + t191 * t142 + t144 * t225 - t147 * t226;
t130 = t144 * t191 - t147 * t188;
t240 = pkin(4) * t246;
t126 = -t130 + t240;
t242 = (t126 + t130 + t240) * qJD(5) - pkin(7) * qJDD(5);
t200 = -g(1) * t149 + g(2) * t148 - t191 * t137 + t188 * t142 + t144 * t226 + t147 * t225;
t181 = qJDD(1) - qJDD(4);
t241 = pkin(4) * t181;
t239 = pkin(1) * qJDD(1);
t207 = t185 * t191 + t186 * t188;
t236 = (t207 * qJD(2) + qJD(4) * t232) * t246;
t235 = (t144 * t188 + t147 * t191) * t246;
t234 = qJDD(3) + g(3);
t230 = t192 * pkin(1) + t189 * qJ(2);
t229 = g(1) * t189 - g(2) * t192;
t183 = t187 ^ 2;
t228 = -t190 ^ 2 + t183;
t223 = qJ(2) * qJDD(1);
t219 = 0.2e1 * t222;
t213 = qJDD(2) - t239;
t210 = g(1) * t192 + g(2) * t189;
t209 = (-t185 * t227 + t165) * t185 - t147 * t186;
t208 = t160 * t191 - t188 * t164;
t156 = t185 * t188 - t186 * t191;
t194 = qJD(5) ^ 2;
t167 = qJDD(5) * t190 - t187 * t194;
t166 = qJDD(5) * t187 + t190 * t194;
t206 = -t241 - t200;
t202 = pkin(7) * t181 + t126 * t246 - t199;
t128 = -t156 * qJD(2) + t208 * qJD(4);
t132 = pkin(4) - t208;
t133 = -pkin(7) + t232;
t201 = -qJDD(5) * t133 + (-t132 * t246 - t126 - t128) * qJD(5);
t198 = pkin(7) * t194 - t206 + t235 + t241;
t197 = -t132 * t181 + t133 * t194 + t206 - t236;
t196 = -t181 * MDP(10) + (-t224 * t246 * t247 - t181 * t183) * MDP(13) + 0.2e1 * (qJD(5) * t228 * t246 - t181 * t187 * t190) * MDP(14) + t166 * MDP(15) + t167 * MDP(16);
t195 = qJD(1) ^ 2;
t176 = t192 * qJ(2);
t159 = t185 * t189 + t186 * t192;
t157 = t185 * t192 - t186 * t189;
t141 = t162 + (-t222 - t223) * t185;
t1 = [qJDD(1) * MDP(1) + t229 * MDP(2) + t210 * MDP(3) + (-qJDD(2) + t229 + 0.2e1 * t239) * MDP(4) + (-t210 + t219 + 0.2e1 * t223) * MDP(5) + (-t213 * pkin(1) - g(1) * (-pkin(1) * t189 + t176) - g(2) * t230 + (t219 + t223) * qJ(2)) * MDP(6) + (0.2e1 * t218 - g(1) * t157 - g(2) * t159 - t162 + (-t163 + t237) * qJDD(1)) * MDP(7) + (-g(1) * t159 + g(2) * t157 + qJDD(1) * t164 + 0.2e1 * t217 + t231) * MDP(8) + (t142 * t164 + t141 * t163 - g(1) * (t193 * t189 + t176) - g(2) * (pkin(2) * t192 + t230) - t209 * qJD(2)) * MDP(9) + (-t208 * t181 + t200 + t236) * MDP(11) + (t128 * t246 + t232 * t181 + t199) * MDP(12) + (t201 * t187 - t197 * t190) * MDP(18) + (t197 * t187 + t201 * t190) * MDP(19) - t196; -qJDD(1) * MDP(4) - t195 * MDP(5) + (-qJ(2) * t195 + t213 - t229) * MDP(6) + (-t185 * t195 - t221) * MDP(7) + (qJDD(1) * t185 - t186 * t195) * MDP(8) + (t209 * qJD(1) + t141 * t186 + t142 * t185 - t229) * MDP(9) + (t244 * t181 - (MDP(18) * qJD(5) * t247 - MDP(12) * t246 + 0.2e1 * t224 * MDP(19)) * t246) * t156 + (t181 * MDP(12) - t166 * MDP(18) - t167 * MDP(19) - t244 * t180) * t207; t167 * MDP(18) - t166 * MDP(19) + t234 * MDP(9); (-t200 - t235) * MDP(11) + (-t130 * t246 - t199) * MDP(12) + (-t198 * MDP(18) + MDP(19) * t242) * t190 + (MDP(18) * t242 + t198 * MDP(19)) * t187 + t196; qJDD(5) * MDP(17) + t228 * MDP(14) * t180 + (-t181 * MDP(16) + t234 * MDP(18) + t202 * MDP(19)) * t190 + (-t180 * t190 * MDP(13) - t181 * MDP(15) + t202 * MDP(18) - t234 * MDP(19)) * t187;];
tau = t1;
