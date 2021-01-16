% Calculate vector of inverse dynamics joint torques for
% S5RPPRP3
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:04:43
% EndTime: 2021-01-15 17:04:46
% DurationCPUTime: 1.11s
% Computational Cost: add. (647->191), mult. (1017->213), div. (0->0), fcn. (491->8), ass. (0->93)
t172 = cos(qJ(4));
t168 = cos(pkin(7));
t150 = -pkin(1) * t168 - pkin(2);
t143 = -pkin(6) + t150;
t133 = qJD(1) * t143 + qJD(3);
t228 = -qJ(5) * qJD(1) + t133;
t231 = t228 * t172;
t163 = qJD(3) * qJD(1);
t161 = qJ(1) + pkin(7);
t155 = sin(t161);
t156 = cos(t161);
t187 = -g(1) * t156 - g(2) * t155;
t167 = sin(pkin(7));
t147 = pkin(1) * t167 + qJ(3);
t214 = qJDD(1) * t147;
t179 = t187 + t214;
t230 = t163 + t179;
t229 = qJDD(2) - g(3);
t149 = g(2) * t156;
t222 = g(1) * t155;
t194 = -t149 + t222;
t170 = sin(qJ(4));
t164 = t170 ^ 2;
t165 = t172 ^ 2;
t227 = MDP(9) * (t164 - t165);
t216 = qJ(5) - t143;
t128 = t216 * t170;
t135 = pkin(4) * t170 + t147;
t211 = qJD(1) * t135;
t127 = qJD(5) + t211;
t226 = -qJD(1) * t127 - t194;
t225 = (qJD(5) + t127) * qJD(1);
t206 = qJD(1) * qJD(4);
t195 = t172 * t206;
t203 = qJDD(1) * t170;
t224 = qJDD(5) + (t195 + t203) * pkin(4);
t171 = sin(qJ(1));
t221 = g(1) * t171;
t159 = g(3) * t170;
t173 = cos(qJ(1));
t158 = t173 * pkin(1);
t220 = qJD(4) * pkin(4);
t219 = qJDD(4) * pkin(4);
t218 = t133 * t170;
t175 = qJD(1) ^ 2;
t217 = t170 * t175;
t120 = -qJD(2) * t170 + t231;
t119 = t120 + t220;
t215 = t119 - t120;
t138 = qJD(1) * t147;
t210 = qJD(1) * t138;
t209 = qJD(1) * t170;
t207 = qJ(5) * qJDD(1);
t205 = qJD(1) * qJD(5);
t204 = qJD(2) * qJD(4);
t202 = qJDD(1) * t172;
t201 = qJDD(4) * t143;
t200 = qJDD(4) * t170;
t199 = MDP(13) + MDP(15);
t198 = MDP(14) + MDP(16);
t197 = t156 * pkin(2) + t155 * qJ(3) + t158;
t134 = t163 + t214;
t196 = t170 * t206;
t129 = t216 * t172;
t193 = (-t164 - t165) * MDP(17);
t191 = t127 + t211;
t190 = 0.2e1 * t138;
t132 = qJDD(1) * t143 + qJDD(3);
t189 = -t132 + t207;
t174 = qJD(4) ^ 2;
t137 = qJDD(4) * t172 - t174 * t170;
t188 = t150 * MDP(7);
t185 = -t210 - t222;
t121 = -qJ(5) * t209 + qJD(2) * t172 + t218;
t184 = t119 * t172 + t121 * t170;
t183 = -pkin(1) * t171 - pkin(2) * t155 + t156 * qJ(3);
t182 = -t204 - t207;
t125 = t172 * t132;
t181 = -t170 * qJDD(2) + t172 * t149 + t125 + t159;
t152 = t170 * t204;
t180 = t170 * t222 - t229 * t172 + t152;
t122 = t134 + t224;
t141 = t172 * t220 + qJD(3);
t178 = qJD(1) * t141 + qJDD(1) * t135 + t122 + t187;
t177 = -t143 * t174 + t134 + t230;
t169 = -qJ(5) - pkin(6);
t142 = qJ(5) * t196;
t136 = -t172 * t174 - t200;
t124 = -qJD(4) * t129 - qJD(5) * t170;
t123 = qJD(4) * t128 - qJD(5) * t172;
t118 = -t152 + (qJD(4) * t228 + qJDD(2)) * t172 + (-t189 - t205) * t170;
t117 = t219 + t125 + t142 + (-qJD(4) * t133 - qJDD(2)) * t170 + (t182 - t205) * t172;
t1 = [(-g(2) * t173 + t221) * MDP(2) + (g(1) * t173 + g(2) * t171) * MDP(3) + (pkin(1) * t221 - g(2) * t158) * MDP(4) + (qJDD(3) - t194) * MDP(5) + (0.2e1 * t163 + t179) * MDP(6) + (-g(1) * t183 - g(2) * t197 + t138 * qJD(3) + qJDD(3) * t150 + t134 * t147) * MDP(7) + t137 * MDP(10) + t136 * MDP(11) + (qJD(4) * t123 - qJDD(4) * t129) * MDP(15) + (-qJD(4) * t124 + qJDD(4) * t128) * MDP(16) + t194 * MDP(17) + (-t118 * t128 + t121 * t124 - t117 * t129 + t119 * t123 + t122 * t135 + t127 * t141 - g(1) * (t155 * t169 + t183) - g(2) * (-t156 * t169 + t197)) * MDP(18) + 0.2e1 * t206 * t227 + (MDP(1) + t147 * MDP(6) + t165 * MDP(8) + (t167 ^ 2 + t168 ^ 2) * MDP(4) * pkin(1) ^ 2 + (0.2e1 * MDP(5) + t188) * t150) * qJDD(1) + (t177 * MDP(13) - MDP(14) * t201 + t178 * MDP(15) + (-qJD(1) * t124 + qJDD(1) * t128 - t118) * MDP(17) + t187 * MDP(18) * pkin(4) + (-t190 * MDP(14) - t191 * MDP(16) + (-qJD(1) * t129 + t119) * MDP(17)) * qJD(4)) * t170 + (-0.2e1 * MDP(9) * t203 + MDP(13) * t201 + t177 * MDP(14) + t178 * MDP(16) + (-qJD(1) * t123 + qJDD(1) * t129 - t117) * MDP(17) + (-0.2e1 * MDP(8) * t209 + t190 * MDP(13) + t191 * MDP(15) + (qJD(1) * t128 - t121) * MDP(17)) * qJD(4)) * t172; (-qJD(4) * t184 - t117 * t170 + t118 * t172 - g(3)) * MDP(18) + (MDP(4) + MDP(7)) * t229 - t198 * t137 + t199 * t136; -t175 * MDP(6) + (qJDD(3) + t149 + t185) * MDP(7) + (t117 * t172 + t118 * t170 + (-t119 * t170 + t121 * t172) * qJD(4) + t226) * MDP(18) + t199 * (t137 - t217) + t198 * (-t200 + (-t174 - t175) * t172) + (MDP(5) + t188 + t193) * qJDD(1); t172 * MDP(8) * t217 - t175 * t227 + MDP(10) * t202 - MDP(11) * t203 + qJDD(4) * MDP(12) + (t172 * t185 + t181) * MDP(13) + ((-t132 - t204 + t210 - t149) * t170 + t180) * MDP(14) + (0.2e1 * t219 + t142 + (t121 - t218) * qJD(4) + (-pkin(4) * t217 + t182 - t222 - t225) * t172 + t181) * MDP(15) + (-pkin(4) * t165 * t175 + (t120 - t231) * qJD(4) + (-t149 + t189 + t225) * t170 + t180) * MDP(16) + (-pkin(4) * t202 + (-t215 + t220) * t209) * MDP(17) + (t215 * t121 + (t226 * t172 + t117 + t159) * pkin(4)) * MDP(18); (0.2e1 * t195 + t203) * MDP(15) + (-0.2e1 * t196 + t202) * MDP(16) + (qJD(1) * t184 + t224 + t230) * MDP(18) + t175 * t193;];
tau = t1;
