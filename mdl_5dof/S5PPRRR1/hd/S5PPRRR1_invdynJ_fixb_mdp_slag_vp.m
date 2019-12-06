% Calculate vector of inverse dynamics joint torques for
% S5PPRRR1
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PPRRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:54
% EndTime: 2019-12-05 15:12:56
% DurationCPUTime: 0.74s
% Computational Cost: add. (613->130), mult. (1176->184), div. (0->0), fcn. (960->14), ass. (0->76)
t173 = sin(pkin(8));
t175 = cos(pkin(8));
t225 = g(1) * t175 + g(2) * t173;
t172 = sin(pkin(9));
t174 = cos(pkin(9));
t178 = sin(qJ(3));
t181 = cos(qJ(3));
t143 = t172 * t181 + t174 * t178;
t141 = t143 * qJD(3);
t205 = qJDD(1) * t181;
t206 = qJDD(1) * t178;
t198 = -t172 * t206 + t174 * t205;
t130 = qJDD(3) * pkin(3) - qJD(1) * t141 + t198;
t228 = -t172 * t178 + t174 * t181;
t140 = t228 * qJD(3);
t131 = qJD(1) * t140 + t143 * qJDD(1);
t168 = pkin(9) + qJ(3);
t165 = qJ(4) + t168;
t159 = cos(t165);
t177 = sin(qJ(4));
t180 = cos(qJ(4));
t232 = -g(3) * t159 + t180 * t130 - t131 * t177;
t167 = qJDD(3) + qJDD(4);
t220 = pkin(4) * t167;
t138 = t228 * qJD(1);
t137 = qJD(3) * pkin(3) + t138;
t139 = t143 * qJD(1);
t213 = t139 * t180;
t126 = t137 * t177 + t213;
t224 = t126 * qJD(4);
t231 = -t220 + t224 - t232;
t158 = sin(t165);
t230 = t225 * t158;
t209 = qJD(4) * t177;
t192 = g(3) * t158 - t130 * t177 + t139 * t209 + t225 * t159;
t229 = t192 - (qJD(4) * t137 + t131) * t180;
t214 = t139 * t177;
t129 = t138 * t180 - t214;
t222 = pkin(3) * t177;
t160 = pkin(7) + t222;
t221 = pkin(3) * t180;
t161 = -pkin(4) - t221;
t169 = qJD(3) + qJD(4);
t226 = -qJDD(5) * t160 + (-qJD(4) * t221 + t161 * t169 + t129) * qJD(5);
t223 = pkin(3) * t169;
t216 = t126 * t169;
t215 = (t138 * t177 + t213) * t169;
t176 = sin(qJ(5));
t170 = t176 ^ 2;
t179 = cos(qJ(5));
t210 = -t179 ^ 2 + t170;
t208 = qJD(5) * t169;
t207 = qJD(5) * t179;
t203 = -t137 - t223;
t133 = t143 * t180 + t177 * t228;
t196 = -t143 * t177 + t180 * t228;
t197 = -(t133 * qJD(4) + t140 * t177 + t141 * t180) * t169 + t196 * t167;
t125 = t137 * t180 - t214;
t194 = g(1) * t173 - g(2) * t175 - qJDD(2);
t182 = qJD(5) ^ 2;
t191 = -pkin(7) * t182 + t216 + t220;
t190 = t133 * t182 - t197;
t189 = t230 + t232;
t187 = -pkin(4) * t208 - pkin(7) * qJDD(5) + qJD(5) * t125;
t119 = t196 * qJD(4) + t140 * t180 - t141 * t177;
t186 = -qJD(5) * t119 - qJDD(5) * t133 - t196 * t208;
t185 = -t160 * t182 - t161 * t167 - t209 * t223 + t215;
t123 = -pkin(4) * t169 - t125;
t184 = -pkin(7) * t167 - t123 * t169 + t229;
t152 = qJDD(5) * t176 + t179 * t182;
t153 = qJDD(5) * t179 - t176 * t182;
t183 = 0.2e1 * (t167 * t176 * t179 - t210 * t208) * MDP(10) + (0.2e1 * t169 * t176 * t207 + t167 * t170) * MDP(9) + t152 * MDP(11) + t153 * MDP(12) + t167 * MDP(6) + (t123 * t207 + t231 * t176) * MDP(15) + (t123 * qJD(5) * t176 + t230 * t179) * MDP(14);
t166 = t169 ^ 2;
t164 = cos(t168);
t163 = sin(t168);
t1 = [(qJDD(1) - g(3)) * MDP(1) + (-g(3) + (t172 ^ 2 + t174 ^ 2) * qJDD(1)) * MDP(2) + (-qJD(3) * t141 + qJDD(3) * t228) * MDP(4) + (-qJD(3) * t140 - qJDD(3) * t143) * MDP(5) + t197 * MDP(7) + (-t119 * t169 - t133 * t167) * MDP(8) + (-t190 * MDP(14) + t186 * MDP(15)) * t179 + (t186 * MDP(14) + t190 * MDP(15)) * t176; t153 * MDP(14) - t152 * MDP(15) - t194 * MDP(2); qJDD(3) * MDP(3) + (-g(3) * t164 + t225 * t163 + t198) * MDP(4) + (g(3) * t163 + t225 * t164 - t172 * t205 - t174 * t206) * MDP(5) + (t167 * t221 + t189 + t215) * MDP(7) + (t129 * t169 - t131 * t180 - t167 * t222 + t192) * MDP(8) + (t203 * MDP(7) * t177 + (-t139 * MDP(7) + t203 * MDP(8)) * t180) * qJD(4) + (t139 * MDP(4) + t138 * MDP(5) + (-t143 * MDP(4) - MDP(5) * t228) * qJD(1)) * qJD(3) + ((t185 - t231) * MDP(14) + t226 * MDP(15)) * t179 + ((-t185 - t230) * MDP(15) + t226 * MDP(14)) * t176 + t183; (t189 + t216 - t224) * MDP(7) + (t125 * t169 + t229) * MDP(8) + ((t191 - t231) * MDP(14) + t187 * MDP(15)) * t179 + (t187 * MDP(14) + (-t191 - t230) * MDP(15)) * t176 + t183; qJDD(5) * MDP(13) + t210 * MDP(10) * t166 + (t167 * MDP(12) - t194 * MDP(14) + t184 * MDP(15)) * t179 + (-t166 * t179 * MDP(9) + t167 * MDP(11) + t184 * MDP(14) + t194 * MDP(15)) * t176;];
tau = t1;
