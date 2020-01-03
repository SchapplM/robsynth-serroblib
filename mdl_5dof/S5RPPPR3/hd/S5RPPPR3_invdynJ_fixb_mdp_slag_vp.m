% Calculate vector of inverse dynamics joint torques for
% S5RPPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:05
% EndTime: 2019-12-31 17:44:07
% DurationCPUTime: 1.00s
% Computational Cost: add. (513->179), mult. (943->238), div. (0->0), fcn. (633->10), ass. (0->92)
t185 = sin(pkin(8));
t182 = t185 ^ 2;
t187 = cos(pkin(8));
t236 = t187 ^ 2 + t182;
t186 = sin(pkin(7));
t171 = pkin(1) * t186 + qJ(3);
t154 = qJD(1) * qJD(3) + qJDD(1) * t171;
t235 = MDP(8) + MDP(12);
t184 = qJ(1) + pkin(7);
t180 = cos(t184);
t175 = g(2) * t180;
t179 = sin(t184);
t232 = g(1) * t179;
t234 = -t175 + t232;
t233 = pkin(3) * t187;
t231 = -pkin(6) + t171;
t134 = t185 * qJDD(2) + t187 * t154;
t230 = t134 * t185 - g(3);
t229 = pkin(6) * qJDD(1);
t128 = t134 * t187;
t163 = t171 * qJD(1);
t228 = (t185 * qJD(2) + t187 * t163) * t187;
t227 = t179 * t187;
t191 = cos(qJ(5));
t226 = t185 * t191;
t189 = sin(qJ(5));
t225 = t187 * t189;
t224 = qJD(3) * t228 + t171 * t128;
t223 = MDP(18) * t191;
t222 = qJD(1) * t189;
t221 = qJD(1) * t191;
t220 = qJD(4) * t185;
t218 = qJD(1) * qJD(4);
t188 = cos(pkin(7));
t176 = -pkin(1) * t188 - pkin(2);
t200 = qJ(4) * t185 - t176;
t151 = -t200 - t233;
t217 = qJDD(1) * t151;
t215 = qJDD(1) * t176;
t214 = qJDD(1) * t189;
t213 = qJDD(1) * t191;
t192 = cos(qJ(1));
t212 = t192 * pkin(1) + t180 * pkin(2) + t179 * qJ(3);
t211 = t185 * t221;
t210 = qJD(5) * t226;
t209 = t187 * t222;
t208 = t185 * t218;
t142 = qJD(2) * t187 - t185 * t163;
t133 = qJDD(2) * t187 - t185 * t154;
t207 = -g(1) * t180 - g(2) * t179;
t190 = sin(qJ(1));
t206 = g(1) * t190 - g(2) * t192;
t166 = t185 * t213;
t204 = -t187 * t214 + t166;
t203 = qJDD(3) - t208;
t152 = t231 * t185;
t153 = t231 * t187;
t202 = t152 * t191 - t153 * t189;
t201 = t152 * t189 + t153 * t191;
t158 = -t225 + t226;
t157 = t185 * t189 + t187 * t191;
t131 = qJDD(4) - t133;
t199 = -pkin(1) * t190 - pkin(2) * t179 + t180 * qJ(3);
t124 = t203 + t217;
t198 = -t124 - t217 - t175;
t161 = qJDD(3) + t215;
t197 = -t161 - t215 - t175;
t146 = t157 * qJD(5);
t141 = (pkin(3) + pkin(4)) * t187 + t200;
t196 = t236 * t154 + t128 + t207;
t194 = qJD(1) ^ 2;
t193 = qJD(5) ^ 2;
t165 = qJD(5) * t209;
t164 = g(1) * t227;
t150 = -t209 + t211;
t148 = t157 * qJD(1);
t147 = -qJD(5) * t225 + t210;
t140 = qJD(4) - t142;
t139 = t157 * t180;
t138 = t158 * t180;
t137 = t157 * t179;
t136 = t158 * t179;
t135 = qJD(1) * t151 + qJD(3);
t126 = qJD(1) * t141 - qJD(3);
t125 = -t187 * t229 + t134;
t123 = -t185 * t229 + t131;
t121 = -qJD(5) * t147 - qJDD(5) * t157;
t120 = -qJD(5) * t146 + qJDD(5) * t158;
t119 = qJD(1) * t210 + qJDD(1) * t157 - t165;
t118 = -qJD(1) * t146 + t204;
t117 = qJDD(1) * t141 - t203;
t1 = [qJDD(1) * MDP(1) + t206 * MDP(2) + (g(1) * t192 + g(2) * t190) * MDP(3) + (t206 + (t186 ^ 2 + t188 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t187 * t197 + t164) * MDP(5) + (-t197 - t232) * t185 * MDP(6) + (-t133 * t185 + t196) * MDP(7) + (t161 * t176 - g(1) * t199 - g(2) * t212 + (-qJD(3) * t142 - t133 * t171) * t185 + t224) * MDP(8) + (t164 + (t198 + t208) * t187) * MDP(9) + (t131 * t185 + t196) * MDP(10) + (t182 * t218 + (t198 + t232) * t185) * MDP(11) + (t124 * t151 - g(1) * (-pkin(3) * t227 + t199) - g(2) * (t180 * t233 + t212) + (t234 * qJ(4) + t140 * qJD(3) - t135 * qJD(4) + t131 * t171) * t185 + t224) * MDP(12) + (t118 * t158 - t146 * t150) * MDP(13) + (-t118 * t157 - t119 * t158 + t146 * t148 - t147 * t150) * MDP(14) + t120 * MDP(15) + t121 * MDP(16) + (t148 * t220 + t141 * t119 + t117 * t157 + t126 * t147 + t202 * qJDD(5) + g(1) * t137 - g(2) * t139 + (qJD(3) * t158 - qJD(5) * t201) * qJD(5)) * MDP(18) + (t150 * t220 + t141 * t118 + t117 * t158 - t126 * t146 - t201 * qJDD(5) + g(1) * t136 - g(2) * t138 + (-qJD(3) * t157 - qJD(5) * t202) * qJD(5)) * MDP(19); (qJDD(2) - g(3)) * MDP(4) + (t133 * t187 + t230) * MDP(8) + (-t131 * t187 + t230) * MDP(12) + t121 * MDP(18) - t120 * MDP(19); (-qJD(5) * t150 + t165) * MDP(18) + (qJD(5) * t148 - t166) * MDP(19) - (MDP(7) + MDP(10)) * t236 * t194 + ((t142 * t185 - t228) * MDP(8) + (-t140 * t185 - t220 - t228) * MDP(12) + (MDP(19) * t157 - t185 * t223) * qJD(5)) * qJD(1) + ((-MDP(12) * qJ(4) - MDP(18) * t189 - MDP(11) + MDP(6)) * t185 + (-MDP(12) * pkin(3) + MDP(19) * t189 - MDP(5) - MDP(9) - t223) * t187 + t235 * t176) * qJDD(1) + t235 * (qJDD(3) - t234); -t182 * t194 * MDP(11) + (g(3) * t187 + t131) * MDP(12) + (qJDD(5) * t191 - t189 * t193) * MDP(18) + (-qJDD(5) * t189 - t191 * t193) * MDP(19) + (-t194 * t187 * MDP(9) + qJDD(1) * MDP(10) + t207 * MDP(12) + (t135 * MDP(12) - t148 * MDP(18) - t150 * MDP(19)) * qJD(1)) * t185; t150 * t148 * MDP(13) + (-t148 ^ 2 + t150 ^ 2) * MDP(14) + t204 * MDP(15) + (-t185 * t214 - t187 * t213 + t165) * MDP(16) + qJDD(5) * MDP(17) + (-g(1) * t138 - g(2) * t136 + g(3) * t157 + t191 * t123 - t189 * t125 - t126 * t150) * MDP(18) + (g(1) * t139 + g(2) * t137 + g(3) * t158 - t189 * t123 - t191 * t125 + t126 * t148) * MDP(19) + ((-t185 * t222 - t187 * t221 + t148) * MDP(15) + (t150 - t211) * MDP(16)) * qJD(5);];
tau = t1;
