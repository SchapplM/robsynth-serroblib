% Calculate vector of inverse dynamics joint torques for
% S4RRRP2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:12
% EndTime: 2019-12-31 17:13:13
% DurationCPUTime: 0.82s
% Computational Cost: add. (605->154), mult. (905->210), div. (0->0), fcn. (439->8), ass. (0->81)
t213 = qJ(4) + pkin(6);
t160 = sin(qJ(2));
t202 = pkin(1) * qJD(1);
t189 = t160 * t202;
t163 = cos(qJ(2));
t207 = pkin(1) * t163;
t199 = qJD(2) * t189 - qJDD(1) * t207;
t153 = qJDD(1) + qJDD(2);
t204 = t153 * pkin(2);
t157 = qJ(1) + qJ(2);
t150 = cos(t157);
t205 = g(2) * t150;
t212 = t199 - t204 + t205;
t192 = qJDD(1) * t160;
t194 = qJD(2) * t163;
t170 = qJD(1) * t194 + t192;
t122 = t170 * pkin(1) + pkin(6) * t153;
t188 = t163 * t202;
t154 = qJD(1) + qJD(2);
t206 = pkin(2) * t154;
t132 = -t188 - t206;
t149 = sin(t157);
t209 = g(1) * t150 + g(2) * t149;
t211 = -t132 * t154 - t122 + t209;
t159 = sin(qJ(3));
t181 = t213 * t154 + t189;
t116 = t181 * t159;
t115 = qJD(3) * pkin(3) - t116;
t141 = g(1) * t149;
t208 = t141 - t205;
t203 = MDP(15) * pkin(3);
t143 = pkin(1) * t160 + pkin(6);
t200 = -qJ(4) - t143;
t155 = t159 ^ 2;
t162 = cos(qJ(3));
t156 = t162 ^ 2;
t197 = t155 - t156;
t144 = pkin(3) * t162 + pkin(2);
t120 = -t144 * t154 + qJD(4) - t188;
t196 = MDP(15) * t120;
t195 = qJD(2) * t160;
t193 = qJD(3) * t159;
t191 = MDP(12) * qJDD(3);
t190 = MDP(13) * qJDD(3);
t187 = pkin(1) * t194;
t186 = pkin(3) * t196;
t185 = t154 * t195;
t184 = t154 * t193;
t182 = qJD(3) * t213;
t180 = qJD(3) * t200;
t179 = (-t155 - t156) * MDP(14);
t178 = MDP(15) * t188;
t161 = sin(qJ(1));
t164 = cos(qJ(1));
t176 = g(1) * t161 - g(2) * t164;
t175 = qJD(3) * t181;
t173 = t188 - t206;
t145 = -pkin(2) - t207;
t172 = t145 * t154 - t187;
t171 = qJ(4) * t153 + qJD(4) * t154 + t122;
t165 = qJD(3) ^ 2;
t169 = pkin(6) * t165 - t154 * t189 - t204;
t111 = pkin(3) * t184 - t144 * t153 + qJDD(4) + t199;
t168 = pkin(1) * t185 + t143 * t165 + t145 * t153;
t167 = -g(1) * (-t144 * t149 + t150 * t213) - g(2) * (t150 * t144 + t149 * t213);
t110 = -t159 * t175 + t171 * t162;
t166 = 0.2e1 * (-t197 * t154 * qJD(3) + t153 * t159 * t162) * MDP(8) + (t153 * t155 + 0.2e1 * t162 * t184) * MDP(7) + (qJDD(3) * t162 - t159 * t165) * MDP(10) + (qJDD(3) * t159 + t162 * t165) * MDP(9) + t153 * MDP(4) + (t132 * t193 + t162 * t141) * MDP(12) + t209 * MDP(6) + (t110 * t162 - t209) * MDP(14) + (t132 * qJD(3) * t162 + t212 * t159) * MDP(13) + (-t199 + t208) * MDP(5);
t152 = t154 ^ 2;
t151 = t162 * qJ(4);
t148 = t162 * qJD(4);
t135 = pkin(6) * t162 + t151;
t134 = t213 * t159;
t128 = t143 * t162 + t151;
t127 = t200 * t159;
t124 = -qJD(4) * t159 - t162 * t182;
t123 = -t159 * t182 + t148;
t117 = t181 * t162;
t114 = (-qJD(4) - t187) * t159 + t162 * t180;
t113 = t159 * t180 + t162 * t187 + t148;
t109 = qJDD(3) * pkin(3) - t171 * t159 - t162 * t175;
t1 = [(-t143 * t191 + (t168 - t141) * MDP(13) + (-t114 * t154 - t127 * t153 - t109) * MDP(14) + (t172 * MDP(12) + (-t128 * t154 - t117) * MDP(14) + t186) * qJD(3)) * t159 + qJDD(1) * MDP(1) + (-t111 * pkin(2) + t109 * t127 + t110 * t128 + t117 * t113 + t115 * t114 + t167) * MDP(15) + ((t153 * t163 - t185) * MDP(5) + (-t153 * t160 - t154 * t194 - t170) * MDP(6) + (-t111 * t163 + t120 * t195 + t176) * MDP(15)) * pkin(1) + t166 + ((-t168 - t212) * MDP(12) - t143 * t190 + (t113 * t154 + t128 * t153) * MDP(14) - t111 * t203 + (t172 * MDP(13) + (-t127 * t154 - t115) * MDP(14)) * qJD(3)) * t162 + t176 * MDP(2) + (g(1) * t164 + g(2) * t161) * MDP(3); (-pkin(6) * t191 + (t169 - t141) * MDP(13) + (-t124 * t154 + t134 * t153 - t109) * MDP(14) + t115 * t178 + (t173 * MDP(12) + (-t135 * t154 - t117) * MDP(14) + t186) * qJD(3)) * t159 + ((-t169 - t212) * MDP(12) - pkin(6) * t190 + (t123 * t154 + t135 * t153) * MDP(14) - t117 * t178 + (t173 * MDP(13) + (t134 * t154 - t115) * MDP(14)) * qJD(3)) * t162 + (-MDP(6) * t192 + ((MDP(5) * t154 - t196) * t160 + (-qJD(2) * MDP(6) + (MDP(6) + t179) * t154) * t163) * qJD(1)) * pkin(1) + (-t109 * t134 + t110 * t135 - t111 * t144 + t115 * t124 + t117 * t123 + t167) * MDP(15) + t166; qJDD(3) * MDP(11) + (pkin(3) * t109 + (t115 + t116) * t117) * MDP(15) + t197 * MDP(8) * t152 + (t153 * MDP(10) + (-MDP(12) - t203) * g(3) + t211 * MDP(13)) * t162 + (-t152 * t162 * MDP(7) + t153 * MDP(9) + t211 * MDP(12) + g(3) * MDP(13) + (-t153 * MDP(14) + (-t120 * t154 + t209) * MDP(15)) * pkin(3)) * t159; ((t115 * t159 - t117 * t162) * t154 + t111 - t208) * MDP(15) + t152 * t179;];
tau = t1;
