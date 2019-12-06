% Calculate vector of inverse dynamics joint torques for
% S5PRRRR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:23
% EndTime: 2019-12-05 17:06:25
% DurationCPUTime: 0.80s
% Computational Cost: add. (712->153), mult. (1008->202), div. (0->0), fcn. (535->12), ass. (0->89)
t185 = cos(qJ(3));
t239 = pkin(2) * t185;
t167 = qJDD(2) * t239;
t175 = qJDD(2) + qJDD(3);
t182 = sin(qJ(3));
t236 = pkin(2) * qJD(2);
t214 = t182 * t236;
t137 = pkin(3) * t175 - qJD(3) * t214 + t167;
t177 = qJD(2) + qJD(3);
t146 = pkin(3) * t177 + t185 * t236;
t176 = pkin(9) + qJ(2);
t174 = qJ(3) + t176;
t163 = qJ(4) + t174;
t158 = cos(t163);
t184 = cos(qJ(4));
t181 = sin(qJ(4));
t221 = qJD(4) * t181;
t244 = -g(2) * t158 + t184 * t137 - t146 * t221;
t219 = qJD(4) * t184;
t211 = t182 * t219;
t216 = qJDD(2) * t182;
t223 = qJD(3) * t185;
t170 = qJDD(4) + t175;
t238 = pkin(4) * t170;
t243 = -t238 + (t181 * t216 + (t181 * t223 + t211) * qJD(2)) * pkin(2) - t244;
t210 = qJD(2) * t223;
t197 = t210 + t216;
t157 = sin(t163);
t205 = t181 * t214;
t199 = g(1) * t158 + g(2) * t157 + qJD(4) * t205 - t181 * t137;
t242 = -(t197 * pkin(2) + qJD(4) * t146) * t184 + t199;
t230 = t181 * t182;
t200 = t184 * t185 - t230;
t141 = t200 * t236;
t164 = pkin(3) * t181 + pkin(8);
t165 = -pkin(3) * t184 - pkin(4);
t173 = qJD(4) + t177;
t241 = -qJDD(5) * t164 + (-pkin(3) * t219 + t165 * t173 + t141) * qJD(5);
t152 = g(1) * t157;
t235 = MDP(9) * t181;
t166 = pkin(3) + t239;
t229 = t182 * t184;
t201 = t181 * t185 + t229;
t234 = (t166 * t221 + (qJD(3) * t201 + t211) * pkin(2)) * t173;
t233 = (t146 * t181 + t184 * t214) * t173;
t232 = t201 * t236 * t173;
t231 = t170 * t184;
t228 = qJDD(1) - g(3);
t134 = t146 * t184 - t205;
t131 = -pkin(4) * t173 - t134;
t180 = sin(qJ(5));
t183 = cos(qJ(5));
t227 = t131 * qJD(5) * t180 + t183 * t152;
t161 = sin(t174);
t162 = cos(t174);
t226 = g(1) * t162 + g(2) * t161;
t225 = pkin(2) * t229 + t181 * t166;
t178 = t180 ^ 2;
t224 = -t183 ^ 2 + t178;
t220 = qJD(4) * t182;
t218 = qJD(5) * t173;
t217 = qJD(5) * t183;
t213 = t131 * t217 + t243 * t180;
t212 = t173 * t221;
t209 = qJD(2) * t220;
t206 = qJD(3) * (-qJD(2) - t177);
t204 = g(1) * t161 - g(2) * t162 + t167;
t186 = qJD(5) ^ 2;
t147 = qJDD(5) * t180 + t183 * t186;
t148 = qJDD(5) * t183 - t180 * t186;
t202 = 0.2e1 * (t170 * t180 * t183 - t218 * t224) * MDP(12) + (0.2e1 * t173 * t180 * t217 + t170 * t178) * MDP(11) + t147 * MDP(13) + t148 * MDP(14) + t170 * MDP(8);
t198 = -MDP(10) * t184 - MDP(7) - t235;
t196 = pkin(8) * t186 - t233 - t238;
t138 = pkin(2) * t230 - t166 * t184 - pkin(4);
t139 = pkin(8) + t225;
t195 = t138 * t170 + t139 * t186 + t234;
t194 = t152 + t244;
t193 = -pkin(8) * t170 - t131 * t173 + t242;
t126 = t166 * t219 + (qJD(3) * t200 - t181 * t220) * pkin(2);
t191 = -qJDD(5) * t139 + (t138 * t173 - t126) * qJD(5);
t190 = -pkin(4) * t218 - pkin(8) * qJDD(5) + qJD(5) * t134;
t189 = -t146 * t219 + t199;
t188 = pkin(3) * t212 + t164 * t186 + t165 * t170 - t232;
t187 = MDP(16) * t227 + MDP(17) * t213 + t202;
t172 = cos(t176);
t171 = sin(t176);
t169 = t173 ^ 2;
t168 = t175 * MDP(5);
t1 = [MDP(1) * t228 + t148 * MDP(16) - t147 * MDP(17); qJDD(2) * MDP(2) + (g(1) * t171 - g(2) * t172) * MDP(3) + (g(1) * t172 + g(2) * t171) * MDP(4) + t168 + ((t175 * t185 + t182 * t206) * pkin(2) + t204) * MDP(6) + (((-qJDD(2) - t175) * t182 + t185 * t206) * pkin(2) + t226) * MDP(7) + (t166 * t231 - t234 + (-t184 * t209 + (-t210 + (-qJDD(2) - t170) * t182) * t181) * pkin(2) + t194) * MDP(9) + (-t126 * t173 - t225 * t170 + t242) * MDP(10) + (t191 * t180 + (-t195 - t243) * t183 + t227) * MDP(16) + (t191 * t183 + (t195 - t152) * t180 + t213) * MDP(17) + t202; t168 + t204 * MDP(6) + t226 * MDP(7) + (t194 + t232) * MDP(9) + (t141 * t173 + t189) * MDP(10) + ((-t212 + t231) * MDP(9) + (-t170 * t181 - t173 * t219) * MDP(10)) * pkin(3) + (t198 * t216 + (((-qJD(3) + t177) * MDP(6) - MDP(9) * t219) * t182 + (t177 * MDP(7) + qJD(3) * t198) * t185) * qJD(2)) * pkin(2) + ((t188 - t152) * MDP(17) + t241 * MDP(16)) * t180 + ((-t188 - t243) * MDP(16) + t241 * MDP(17)) * t183 + t187; (t194 + t233) * MDP(9) + (t134 * t173 + t189) * MDP(10) + (-t197 * t235 + (-MDP(10) * t197 - MDP(9) * t209) * t184) * pkin(2) + (t190 * MDP(16) + (t196 - t152) * MDP(17)) * t180 + ((-t196 - t243) * MDP(16) + t190 * MDP(17)) * t183 + t187; qJDD(5) * MDP(15) + t224 * MDP(12) * t169 + (t170 * MDP(14) + MDP(16) * t228 + MDP(17) * t193) * t183 + (-t169 * t183 * MDP(11) + t170 * MDP(13) + MDP(16) * t193 - MDP(17) * t228) * t180;];
tau = t1;
