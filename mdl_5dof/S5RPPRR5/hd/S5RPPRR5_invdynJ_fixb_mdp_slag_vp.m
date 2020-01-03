% Calculate vector of inverse dynamics joint torques for
% S5RPPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPPRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:41
% EndTime: 2019-12-31 17:56:43
% DurationCPUTime: 0.62s
% Computational Cost: add. (608->125), mult. (847->164), div. (0->0), fcn. (455->10), ass. (0->62)
t138 = sin(pkin(8));
t126 = pkin(1) * t138 + qJ(3);
t175 = t126 * qJDD(1);
t141 = sin(qJ(4));
t146 = qJD(5) ^ 2;
t132 = qJDD(1) - qJDD(4);
t144 = cos(qJ(4));
t179 = t132 * t144;
t168 = qJD(1) - qJD(4);
t188 = t168 ^ 2;
t189 = (t146 + t188) * t141 + t179;
t171 = qJD(5) * t168;
t139 = cos(pkin(8));
t127 = -pkin(1) * t139 - pkin(2);
t125 = -pkin(3) + t127;
t176 = t141 * t125 + t144 * t126;
t169 = qJ(1) + pkin(8);
t129 = cos(t169);
t167 = sin(t169);
t109 = -t129 * t144 - t167 * t141;
t110 = t129 * t141 - t167 * t144;
t113 = t125 * qJDD(1) + qJDD(3);
t115 = t125 * qJD(1) + qJD(3);
t170 = qJD(3) * qJD(1);
t116 = t170 + t175;
t120 = t126 * qJD(1);
t172 = qJD(4) * t144;
t173 = qJD(4) * t141;
t151 = g(1) * t109 + g(2) * t110 + t141 * t113 + t115 * t172 + t144 * t116 - t120 * t173;
t183 = pkin(4) * t168;
t98 = t115 * t144 - t120 * t141;
t96 = -t98 + t183;
t186 = (t96 + t98 + t183) * qJD(5) - pkin(7) * qJDD(5);
t152 = -g(1) * t110 + g(2) * t109 - t144 * t113 + t115 * t173 + t141 * t116 + t120 * t172;
t162 = t127 * qJDD(1);
t184 = pkin(4) * t132;
t182 = t168 * (t115 * t141 + t120 * t144);
t180 = (qJD(3) * t141 + qJD(4) * t176) * t168;
t140 = sin(qJ(5));
t143 = cos(qJ(5));
t178 = t140 * t143;
t177 = qJDD(2) - g(3);
t135 = t140 ^ 2;
t174 = -t143 ^ 2 + t135;
t142 = sin(qJ(1));
t145 = cos(qJ(1));
t163 = g(1) * t142 - g(2) * t145;
t161 = t125 * t144 - t126 * t141;
t160 = -t184 - t152;
t156 = pkin(7) * t132 + t168 * t96 - t151;
t155 = g(1) * t167 - g(2) * t129 - qJDD(3);
t154 = -qJDD(5) * t141 + 0.2e1 * t144 * t171;
t100 = qJD(3) * t144 + t161 * qJD(4);
t103 = pkin(4) - t161;
t104 = -pkin(7) + t176;
t153 = -qJDD(5) * t104 + (-t103 * t168 - t100 - t96) * qJD(5);
t150 = pkin(7) * t146 - t160 + t182 + t184;
t149 = -t103 * t132 + t104 * t146 + t160 - t180;
t117 = qJDD(5) * t140 + t143 * t146;
t118 = qJDD(5) * t143 - t146 * t140;
t148 = (-t132 * t135 - 0.2e1 * t171 * t178) * MDP(11) + 0.2e1 * (-t132 * t178 + t174 * t171) * MDP(12) + t117 * MDP(13) + t118 * MDP(14) - t132 * MDP(8);
t1 = [qJDD(1) * MDP(1) + t163 * MDP(2) + (g(1) * t145 + g(2) * t142) * MDP(3) + (t163 + (t138 ^ 2 + t139 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t155 - 0.2e1 * t162) * MDP(5) + (-g(1) * t129 - g(2) * t167 + 0.2e1 * t170 + 0.2e1 * t175) * MDP(6) + (t116 * t126 + t120 * qJD(3) + (qJDD(3) + t162) * t127 - g(1) * (-t142 * pkin(1) - t167 * pkin(2) + t129 * qJ(3)) - g(2) * (t145 * pkin(1) + t129 * pkin(2) + t167 * qJ(3))) * MDP(7) + (-t161 * t132 + t152 + t180) * MDP(9) + (t100 * t168 + t176 * t132 + t151) * MDP(10) + (t153 * t140 - t149 * t143) * MDP(16) + (t149 * t140 + t153 * t143) * MDP(17) - t148; -MDP(16) * t118 + MDP(17) * t117 + (MDP(4) + MDP(7)) * t177; -qJDD(1) * MDP(5) - qJD(1) ^ 2 * MDP(6) + (-t120 * qJD(1) - t155 + t162) * MDP(7) + (-t141 * t188 - t179) * MDP(9) + (t141 * t132 - t144 * t188) * MDP(10) + (t154 * t140 - t189 * t143) * MDP(16) + (t189 * t140 + t154 * t143) * MDP(17); (-t152 - t182) * MDP(9) + (-t168 * t98 - t151) * MDP(10) + (-t150 * MDP(16) + MDP(17) * t186) * t143 + (MDP(16) * t186 + t150 * MDP(17)) * t140 + t148; qJDD(5) * MDP(15) + t174 * MDP(12) * t188 + (-t132 * MDP(14) - t177 * MDP(16) + t156 * MDP(17)) * t143 + (-MDP(11) * t143 * t188 - t132 * MDP(13) + t156 * MDP(16) + t177 * MDP(17)) * t140;];
tau = t1;
