% Calculate vector of inverse dynamics joint torques for
% S5PPRPR1
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:26
% EndTime: 2019-12-05 15:01:31
% DurationCPUTime: 1.40s
% Computational Cost: add. (557->161), mult. (1201->219), div. (0->0), fcn. (968->14), ass. (0->87)
t178 = sin(pkin(9));
t181 = cos(pkin(9));
t220 = t181 * MDP(6);
t246 = -MDP(7) * t178 + t220;
t179 = sin(pkin(8));
t182 = cos(pkin(8));
t185 = sin(qJ(3));
t187 = cos(qJ(3));
t242 = -t179 * t185 + t182 * t187;
t245 = t242 * qJD(1);
t151 = t179 * t187 + t182 * t185;
t143 = t151 * qJD(1);
t177 = pkin(8) + qJ(3);
t170 = sin(t177);
t172 = cos(t177);
t180 = sin(pkin(7));
t183 = cos(pkin(7));
t202 = g(1) * t183 + g(2) * t180;
t191 = -g(3) * t172 + t170 * t202;
t189 = qJD(3) * t143 + t191;
t147 = t151 * qJD(3);
t138 = qJD(3) * qJ(4) + t143;
t133 = qJD(2) * t181 - t138 * t178;
t134 = qJD(2) * t178 + t138 * t181;
t240 = -t133 * t178 + t134 * t181;
t217 = t178 ^ 2 + t181 ^ 2;
t238 = qJD(3) * t217;
t184 = sin(qJ(5));
t186 = cos(qJ(5));
t150 = t178 * t186 + t181 * t184;
t145 = t150 * qJD(5);
t195 = qJD(4) - t245;
t237 = t151 * qJDD(1);
t232 = g(3) * t170;
t236 = t202 * t172 + t232;
t235 = qJD(5) ^ 2;
t230 = pkin(6) + qJ(4);
t228 = pkin(6) * qJDD(3);
t227 = qJDD(3) * pkin(3);
t224 = t172 * t180;
t223 = t172 * t183;
t125 = qJDD(3) * qJ(4) + t237 + (qJD(4) + t245) * qJD(3);
t122 = qJDD(2) * t178 + t125 * t181;
t146 = t242 * qJD(3);
t215 = qJD(3) * t146;
t214 = qJD(3) * t184;
t213 = qJD(3) * t186;
t211 = qJDD(3) * t184;
t210 = qJDD(3) * t186;
t148 = t178 * t184 - t181 * t186;
t209 = qJDD(5) * t148;
t208 = qJDD(5) * t150;
t205 = t181 * t213;
t207 = qJD(5) * t205 + t178 * t210 + t181 * t211;
t164 = -pkin(4) * t181 - pkin(3);
t206 = t178 * t214;
t204 = -g(1) * t180 + g(2) * t183;
t203 = t217 * qJDD(3);
t160 = t181 * t210;
t200 = -t178 * t211 + t160;
t166 = t181 * qJDD(2);
t121 = -t125 * t178 + t166;
t199 = -t121 * t178 + t122 * t181;
t152 = t230 * t178;
t153 = t230 * t181;
t198 = -t152 * t186 - t153 * t184;
t197 = -t152 * t184 + t153 * t186;
t144 = t148 * qJD(5);
t194 = -qJD(1) * t147 + qJDD(1) * t242;
t193 = -t178 * t213 - t181 * t214;
t192 = qJDD(4) - t194;
t176 = pkin(9) + qJ(5);
t171 = cos(t176);
t169 = sin(t176);
t142 = t150 * qJD(3);
t139 = -t205 + t206;
t137 = -qJD(3) * pkin(3) + t195;
t135 = qJD(3) * t164 + t195;
t132 = -qJD(5) * t145 - t209;
t131 = -qJD(5) * t144 + t208;
t130 = qJD(3) * t145 - t200;
t129 = -qJD(5) * t206 + t207;
t126 = t192 - t227;
t123 = qJDD(3) * t164 + t192;
t120 = t181 * t228 + t122;
t119 = t166 + (-t125 - t228) * t178;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (-g(3) + (t179 ^ 2 + t182 ^ 2) * qJDD(1)) * MDP(2) + (-qJDD(3) * t151 - t215) * MDP(5) + (t151 * t203 + t215 * t217) * MDP(8) + (-t126 * t242 + t137 * t147 + t146 * t240 + t199 * t151 - g(3)) * MDP(9) + (-t242 * t130 + t147 * t139 - t146 * t145 + (t148 * t235 - t208) * t151) * MDP(15) + (-t242 * t129 + t147 * t142 + t146 * t144 + (t150 * t235 + t209) * t151) * MDP(16) + (-MDP(4) - t246) * (qJD(3) * t147 - qJDD(3) * t242); (qJDD(2) + t204) * MDP(2) + (t121 * t181 + t122 * t178 + t204) * MDP(9) + t132 * MDP(15) - t131 * MDP(16); qJDD(3) * MDP(3) + (t189 + t194) * MDP(4) + (-t237 + t236) * MDP(5) + (qJ(4) * t203 + t195 * t238 + t199 - t236) * MDP(8) + (-t126 * pkin(3) - t137 * t143 - g(3) * (pkin(3) * t172 + qJ(4) * t170) + (qJ(4) * t122 + t134 * t195) * t181 + (-qJ(4) * t121 - t133 * t195) * t178 + t202 * (pkin(3) * t170 - qJ(4) * t172)) * MDP(9) + (t129 * t150 - t142 * t144) * MDP(10) + (-t129 * t148 - t130 * t150 + t139 * t144 - t142 * t145) * MDP(11) + t131 * MDP(12) + t132 * MDP(13) + (t198 * qJDD(5) + t164 * t130 + t123 * t148 + t135 * t145 - t143 * t139 + t191 * t171 + (-qJD(5) * t197 - t150 * t195) * qJD(5)) * MDP(15) + (-t197 * qJDD(5) + t164 * t129 + t123 * t150 - t135 * t144 - t143 * t142 - t191 * t169 + (-qJD(5) * t198 + t148 * t195) * qJD(5)) * MDP(16) + t246 * (-t126 + t189 + t227); (t192 - t191) * MDP(9) - t160 * MDP(15) + t207 * MDP(16) + (-t220 - pkin(3) * MDP(9) + (MDP(15) * t184 + MDP(7)) * t178) * qJDD(3) + ((t142 - t193) * MDP(15) + (-t139 - t206) * MDP(16)) * qJD(5) + (-MDP(8) * t238 - MDP(9) * t240) * qJD(3); t142 * t139 * MDP(10) + (-t139 ^ 2 + t142 ^ 2) * MDP(11) + t207 * MDP(12) + t200 * MDP(13) + qJDD(5) * MDP(14) + (-t184 * t120 + t186 * t119 - t135 * t142 - g(1) * (-t169 * t223 + t171 * t180) - g(2) * (-t169 * t224 - t171 * t183) + t169 * t232) * MDP(15) + (-t186 * t120 - t184 * t119 + t135 * t139 - g(1) * (-t169 * t180 - t171 * t223) - g(2) * (t169 * t183 - t171 * t224) + t171 * t232) * MDP(16) + ((t139 - t206) * MDP(12) + (t142 + t193) * MDP(13)) * qJD(5);];
tau = t1;
