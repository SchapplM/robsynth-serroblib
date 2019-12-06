% Calculate vector of inverse dynamics joint torques for
% S5PRPRR2
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PRPRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:17
% EndTime: 2019-12-05 15:45:20
% DurationCPUTime: 0.74s
% Computational Cost: add. (736->135), mult. (1339->191), div. (0->0), fcn. (1027->12), ass. (0->80)
t181 = qJ(2) + pkin(9) + qJ(4);
t176 = sin(t181);
t188 = sin(pkin(8));
t190 = cos(pkin(8));
t215 = g(1) * t190 + g(2) * t188;
t242 = t215 * t176;
t196 = cos(qJ(2));
t173 = qJD(2) * pkin(2) + qJD(1) * t196;
t187 = sin(pkin(9));
t189 = cos(pkin(9));
t193 = sin(qJ(2));
t226 = qJD(1) * t193;
t146 = t189 * t173 - t187 * t226;
t144 = qJD(2) * pkin(3) + t146;
t147 = t173 * t187 + t189 * t226;
t192 = sin(qJ(4));
t195 = cos(qJ(4));
t133 = t144 * t192 + t147 * t195;
t180 = t196 * qJDD(1);
t223 = qJD(1) * qJD(2);
t157 = qJDD(2) * pkin(2) - t193 * t223 + t180;
t208 = qJDD(1) * t193 + t196 * t223;
t138 = t189 * t157 - t208 * t187;
t137 = qJDD(2) * pkin(3) + t138;
t139 = t157 * t187 + t208 * t189;
t177 = cos(t181);
t244 = -g(3) * t177 - t133 * qJD(4) + t195 * t137 - t192 * t139;
t183 = qJDD(2) + qJDD(4);
t237 = pkin(4) * t183;
t243 = -t237 - t244;
t178 = pkin(2) * t189 + pkin(3);
t238 = pkin(2) * t187;
t210 = t195 * t178 - t192 * t238;
t151 = -pkin(4) - t210;
t228 = t192 * t178 + t195 * t238;
t152 = pkin(7) + t228;
t197 = qJD(5) ^ 2;
t159 = t187 * t196 + t189 * t193;
t154 = t159 * qJD(1);
t158 = -t187 * t193 + t189 * t196;
t156 = t158 * qJD(1);
t184 = qJD(2) + qJD(4);
t219 = (-t228 * qJD(4) + t154 * t195 + t156 * t192) * t184;
t241 = -t151 * t183 - t152 * t197 + t219;
t232 = t147 * t192;
t201 = g(3) * t176 - (qJD(4) * t144 + t139) * t195 + qJD(4) * t232 - t192 * t137 + t215 * t177;
t233 = t133 * t184;
t231 = qJDD(1) - g(3);
t230 = -t210 * qJD(4) - t154 * t192 + t156 * t195;
t191 = sin(qJ(5));
t185 = t191 ^ 2;
t194 = cos(qJ(5));
t227 = -t194 ^ 2 + t185;
t225 = qJD(5) * t184;
t224 = qJD(5) * t194;
t132 = t195 * t144 - t232;
t130 = -pkin(4) * t184 - t132;
t222 = t130 * t224 + t243 * t191;
t221 = t130 * qJD(5) * t191 + t242 * t194;
t169 = qJDD(5) * t191 + t194 * t197;
t170 = qJDD(5) * t194 - t191 * t197;
t214 = 0.2e1 * (t183 * t191 * t194 - t227 * t225) * MDP(10) + (0.2e1 * t184 * t191 * t224 + t183 * t185) * MDP(9) + t169 * MDP(11) + t170 * MDP(12) + t183 * MDP(6);
t141 = t158 * t192 + t159 * t195;
t153 = t159 * qJD(2);
t155 = t158 * qJD(2);
t211 = t158 * t195 - t159 * t192;
t213 = -(t141 * qJD(4) + t153 * t195 + t155 * t192) * t184 + t211 * t183;
t209 = g(1) * t188 - g(2) * t190 - qJDD(3);
t207 = -pkin(7) * t197 + t233 + t237;
t206 = t141 * t197 - t213;
t205 = -pkin(4) * t225 - pkin(7) * qJDD(5) + qJD(5) * t132;
t204 = -g(3) * t196 + t215 * t193;
t126 = t211 * qJD(4) - t153 * t192 + t155 * t195;
t203 = -qJD(5) * t126 - qJDD(5) * t141 - t211 * t225;
t202 = -qJDD(5) * t152 + (t151 * t184 + t230) * qJD(5);
t200 = -pkin(7) * t183 - t130 * t184 + t201;
t199 = t242 + t244;
t198 = qJD(2) ^ 2;
t182 = t184 ^ 2;
t1 = [t231 * MDP(1) + (qJDD(2) * t196 - t193 * t198) * MDP(3) + (-qJDD(2) * t193 - t196 * t198) * MDP(4) + (t138 * t158 + t139 * t159 - t146 * t153 + t147 * t155 - g(3)) * MDP(5) + t213 * MDP(7) + (-t126 * t184 - t141 * t183) * MDP(8) + (-t206 * MDP(14) + t203 * MDP(15)) * t194 + (t203 * MDP(14) + t206 * MDP(15)) * t191; qJDD(2) * MDP(2) + (t180 + t204) * MDP(3) + (-t231 * t193 + t215 * t196) * MDP(4) + (t146 * t154 - t147 * t156 + (t138 * t189 + t139 * t187 + t204) * pkin(2)) * MDP(5) + (t210 * t183 + t199 + t219) * MDP(7) + (-t228 * t183 + t230 * t184 + t201) * MDP(8) + (t202 * t191 + (-t243 + t241) * t194 + t221) * MDP(14) + (t202 * t194 + (-t241 - t242) * t191 + t222) * MDP(15) + t214; t170 * MDP(14) - t169 * MDP(15) - t209 * MDP(5); (t199 + t233) * MDP(7) + (t132 * t184 + t201) * MDP(8) + t221 * MDP(14) + t222 * MDP(15) + ((t207 - t243) * MDP(14) + t205 * MDP(15)) * t194 + (t205 * MDP(14) + (-t207 - t242) * MDP(15)) * t191 + t214; qJDD(5) * MDP(13) + t227 * MDP(10) * t182 + (t183 * MDP(12) - t209 * MDP(14) + t200 * MDP(15)) * t194 + (-t182 * t194 * MDP(9) + t183 * MDP(11) + t200 * MDP(14) + t209 * MDP(15)) * t191;];
tau = t1;
