% Calculate vector of inverse dynamics joint torques for
% S4PRRP5
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:36
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:36:11
% EndTime: 2021-01-14 22:36:13
% DurationCPUTime: 0.97s
% Computational Cost: add. (406->162), mult. (890->215), div. (0->0), fcn. (518->6), ass. (0->83)
t139 = qJ(4) + pkin(5);
t143 = cos(qJ(2));
t137 = sin(pkin(6));
t138 = cos(pkin(6));
t160 = g(1) * t138 + g(2) * t137;
t152 = t160 * t143;
t141 = sin(qJ(2));
t198 = g(3) * t141;
t206 = t152 + t198;
t140 = sin(qJ(3));
t166 = t141 * qJD(1) + t139 * qJD(2);
t109 = t166 * t140;
t142 = cos(qJ(3));
t110 = t166 * t142;
t176 = qJD(1) * qJD(2);
t119 = qJDD(2) * pkin(5) + qJDD(1) * t141 + t143 * t176;
t163 = -qJ(4) * qJDD(2) - t119;
t150 = qJD(2) * qJD(4) - t163;
t156 = qJD(3) * t166;
t104 = -t140 * t156 + t150 * t142;
t195 = qJD(3) * pkin(3);
t108 = -t109 + t195;
t205 = -qJD(3) * t108 + t104;
t135 = t140 ^ 2;
t136 = t142 ^ 2;
t183 = t135 + t136;
t204 = t183 * MDP(14);
t187 = qJDD(1) - g(3);
t203 = t160 * t141 + t187 * t143;
t202 = pkin(3) * t135;
t201 = pkin(3) * t142;
t197 = g(3) * t143;
t196 = qJD(2) * pkin(2);
t193 = qJDD(3) * pkin(3);
t192 = t108 * t140;
t191 = t140 * t143;
t190 = t141 * t142;
t189 = t142 * t143;
t145 = qJD(2) ^ 2;
t188 = t142 * t145;
t186 = t109 + t108;
t181 = qJD(1) * t143;
t169 = qJD(3) * t181;
t185 = g(3) * t191 + t142 * t169;
t184 = t135 - t136;
t144 = qJD(3) ^ 2;
t182 = t144 + t145;
t180 = qJD(2) * t142;
t178 = qJD(3) * t140;
t134 = pkin(2) + t201;
t117 = -t134 * qJD(2) + qJD(4) - t181;
t177 = t117 * qJD(2);
t175 = qJD(2) * qJD(3);
t174 = qJDD(1) * t143;
t173 = qJDD(2) * t134;
t172 = qJDD(2) * t140;
t171 = qJDD(2) * t142;
t170 = qJDD(3) * t140;
t168 = t140 * t175;
t133 = t141 * t176;
t167 = qJD(3) * t139;
t128 = -t181 - t196;
t165 = -qJD(2) * t128 - t119;
t164 = 0.2e1 * t142 * t175;
t147 = pkin(3) * t168 + qJDD(4) + t133 - t173;
t107 = t147 - t174;
t162 = t107 - t173;
t161 = t142 * t133 + t140 * t169 + t160 * t190;
t159 = -t110 * t142 + t192;
t157 = -MDP(4) + t204;
t155 = -g(1) * (t137 * t142 - t138 * t191) - g(2) * (-t137 * t191 - t138 * t142) + t140 * t198;
t154 = -g(1) * (-t137 * t140 - t138 * t189) - g(2) * (-t137 * t189 + t138 * t140) + g(3) * t190;
t153 = -0.2e1 * qJDD(2) * pkin(2) + pkin(5) * t144 + t133 - t174;
t151 = 0.2e1 * t168 - t171;
t149 = -pkin(5) * qJDD(3) + (t128 - t196) * qJD(3);
t148 = (-qJD(4) - t117) * qJD(2) + t163;
t146 = (-t160 - t176) * t141;
t123 = t139 * t142;
t122 = t139 * t140;
t112 = -qJD(4) * t140 - t142 * t167;
t111 = qJD(4) * t142 - t140 * t167;
t103 = -t150 * t140 - t142 * t156 + t193;
t1 = [t187 * MDP(1) - g(3) * MDP(15) + (MDP(10) + MDP(12)) * (-t151 * t143 + (-t182 * t142 - t170) * t141) + (MDP(11) + MDP(13)) * ((-qJDD(3) * t141 - 0.2e1 * t143 * t175) * t142 + (-qJDD(2) * t143 + t182 * t141) * t140) + (qJDD(2) * MDP(3) + (-qJD(2) * t192 + t110 * t180 - t107) * MDP(15) + t157 * t145) * t143 + (-t145 * MDP(3) + (-t103 * t140 - t110 * t178 + t205 * t142 + t177) * MDP(15) + t157 * qJDD(2)) * t141; qJDD(2) * MDP(2) + t203 * MDP(3) + (-t187 * t141 + t152) * MDP(4) + (qJDD(2) * t135 + t140 * t164) * MDP(5) + 0.2e1 * (t140 * t171 - t184 * t175) * MDP(6) + (t142 * t144 + t170) * MDP(7) + (qJDD(3) * t142 - t140 * t144) * MDP(8) + (t149 * t140 + (-t153 - t197) * t142 + t161) * MDP(10) + (t149 * t142 + (t146 + t153) * t140 + t185) * MDP(11) + (-qJDD(3) * t122 + (-t162 - t197) * t142 + (t112 + (t117 + (-t134 - t201) * qJD(2)) * t140) * qJD(3) + t161) * MDP(12) + (-qJDD(3) * t123 + (t117 * t142 - t111 + (-t134 * t142 + t202) * qJD(2)) * qJD(3) + (t146 + t162) * t140 + t185) * MDP(13) + ((qJDD(2) * t123 + t205) * t142 + (-qJD(3) * t110 + qJDD(2) * t122 - t103) * t140 + (t111 * t142 - t112 * t140 + (t122 * t142 - t123 * t140) * qJD(3) - t183 * t181) * qJD(2) - t206) * MDP(14) + (t104 * t123 + t110 * t111 - t103 * t122 + t108 * t112 - t107 * t134 + t117 * pkin(3) * t178 - g(3) * (t134 * t143 + t139 * t141) + (-t117 * t141 + t159 * t143) * qJD(1) + t160 * (t134 * t141 - t139 * t143)) * MDP(15); -t140 * MDP(5) * t188 + t184 * t145 * MDP(6) + MDP(7) * t172 + MDP(8) * t171 + qJDD(3) * MDP(9) + (t165 * t140 + t155) * MDP(10) + (t165 * t142 + t154) * MDP(11) + (0.2e1 * t193 + (pkin(3) * t188 + t148) * t140 + t155) * MDP(12) + (t148 * t142 - t145 * t202 + t154) * MDP(13) + (-pkin(3) * t172 + (t186 - t195) * t180) * MDP(14) + (t186 * t110 + (t103 + (-g(1) * t137 + g(2) * t138) * t142 + (-t177 + t206) * t140) * pkin(3)) * MDP(15); t151 * MDP(12) + (t164 + t172) * MDP(13) + (t159 * qJD(2) + t147 - t203) * MDP(15) - t145 * t204;];
tau = t1;
