% Calculate vector of inverse dynamics joint torques for
% S4RRPR3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4RRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:36
% EndTime: 2019-12-31 17:01:37
% DurationCPUTime: 0.44s
% Computational Cost: add. (374->105), mult. (654->156), div. (0->0), fcn. (360->12), ass. (0->67)
t148 = qJ(1) + qJ(2);
t140 = pkin(7) + t148;
t131 = cos(t140);
t144 = qJDD(1) + qJDD(2);
t155 = cos(qJ(2));
t186 = t155 * pkin(1);
t138 = qJDD(1) * t186;
t152 = sin(qJ(2));
t170 = pkin(1) * qJD(1) * t152;
t110 = t144 * pkin(2) - qJD(2) * t170 + t138;
t149 = sin(pkin(7));
t150 = cos(pkin(7));
t176 = qJD(1) * t155;
t167 = qJD(2) * t176;
t171 = qJDD(1) * t152;
t99 = -(t167 + t171) * pkin(1) * t149 + t150 * t110;
t193 = -t144 * pkin(3) + g(2) * t131 - t99;
t141 = sin(t148);
t142 = cos(t148);
t192 = g(1) * t141 - g(2) * t142;
t191 = t152 * MDP(5) + t155 * MDP(6);
t130 = sin(t140);
t189 = g(1) * t130;
t185 = t149 * t152;
t137 = pkin(2) + t186;
t183 = t150 * t137;
t182 = t150 * t152;
t179 = qJDD(3) - g(3);
t129 = pkin(1) * t182;
t178 = t149 * t137 + t129;
t151 = sin(qJ(4));
t146 = t151 ^ 2;
t154 = cos(qJ(4));
t177 = -t154 ^ 2 + t146;
t145 = qJD(1) + qJD(2);
t175 = qJD(4) * t145;
t174 = t154 * qJD(4);
t169 = pkin(1) * t176;
t100 = t150 * pkin(1) * t167 + qJDD(1) * t129 + t149 * t110;
t153 = sin(qJ(1));
t156 = cos(qJ(1));
t166 = g(1) * t153 - g(2) * t156;
t165 = pkin(1) * (t150 * t155 - t185);
t120 = t145 * pkin(2) + t169;
t106 = t150 * t120 - t149 * t170;
t111 = pkin(1) * t185 - pkin(3) - t183;
t112 = pkin(6) + t178;
t114 = (t149 * t155 + t182) * qJD(2) * pkin(1);
t157 = qJD(4) ^ 2;
t163 = t111 * t144 + t112 * t157 + t114 * t145;
t128 = t150 * t170;
t113 = t149 * t169 + t128;
t132 = t149 * pkin(2) + pkin(6);
t133 = -t150 * pkin(2) - pkin(3);
t162 = -t113 * t145 + t132 * t157 + t133 * t144;
t104 = -t145 * pkin(3) - t106;
t161 = -t144 * pkin(6) + g(1) * t131 + g(2) * t130 - t104 * t145 - t100;
t116 = qJD(2) * t165;
t160 = -qJD(4) * t116 - qJDD(4) * t112 + t111 * t175;
t115 = qJD(1) * t165;
t159 = qJD(4) * t115 - qJDD(4) * t132 + t133 * t175;
t121 = qJDD(4) * t151 + t157 * t154;
t122 = qJDD(4) * t154 - t157 * t151;
t158 = 0.2e1 * (t151 * t144 * t154 - t175 * t177) * MDP(9) + (0.2e1 * t145 * t151 * t174 + t146 * t144) * MDP(8) + t121 * MDP(10) + t122 * MDP(11) + t144 * MDP(4) + (t104 * qJD(4) * t151 + t154 * t189) * MDP(13) + (g(1) * t142 + g(2) * t141) * MDP(6) + (t104 * t174 + t193 * t151) * MDP(14) + (t138 + t192) * MDP(5);
t143 = t145 ^ 2;
t107 = t149 * t120 + t128;
t1 = [t158 + (t144 * t155 * MDP(5) + t166 * MDP(7) + ((-qJDD(1) - t144) * MDP(6) - t149 * t99 * MDP(7)) * t152 + t191 * qJD(2) * (-qJD(1) - t145)) * pkin(1) + ((-t163 - t193) * MDP(13) + t160 * MDP(14)) * t154 + (t160 * MDP(13) + (t163 - t189) * MDP(14)) * t151 + (t192 * pkin(2) + t100 * t178 - t106 * t114 + t107 * t116 + t183 * t99) * MDP(7) + t166 * MDP(2) + (g(1) * t156 + g(2) * t153) * MDP(3) + qJDD(1) * MDP(1); (-MDP(6) * t171 + t191 * qJD(1) * (-qJD(2) + t145)) * pkin(1) + (t159 * MDP(13) + (t162 - t189) * MDP(14)) * t151 + ((-t162 - t193) * MDP(13) + t159 * MDP(14)) * t154 + t158 + (t106 * t113 - t107 * t115 + (t100 * t149 + t150 * t99 + t192) * pkin(2)) * MDP(7); t122 * MDP(13) - t121 * MDP(14) + MDP(7) * t179; qJDD(4) * MDP(12) + t177 * MDP(9) * t143 + (t144 * MDP(11) + MDP(13) * t179 + MDP(14) * t161) * t154 + (-t143 * t154 * MDP(8) + t144 * MDP(10) + MDP(13) * t161 - MDP(14) * t179) * t151;];
tau = t1;
