% Calculate vector of inverse dynamics joint torques for
% S5PPRRR5
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PPRRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:47
% EndTime: 2019-12-31 17:35:49
% DurationCPUTime: 0.54s
% Computational Cost: add. (434->120), mult. (716->168), div. (0->0), fcn. (514->10), ass. (0->66)
t139 = cos(pkin(8));
t174 = qJ(3) + qJ(4);
t163 = sin(t174);
t164 = cos(t174);
t177 = sin(pkin(8));
t111 = -t139 * t164 - t177 * t163;
t112 = t139 * t163 - t177 * t164;
t161 = g(1) * t112 - g(2) * t111;
t145 = cos(qJ(3));
t132 = t145 * qJDD(2);
t142 = sin(qJ(3));
t168 = qJD(2) * qJD(3);
t115 = qJDD(3) * pkin(3) - t142 * t168 + t132;
t144 = cos(qJ(4));
t110 = t144 * t115;
t141 = sin(qJ(4));
t125 = qJD(3) * pkin(3) + t145 * qJD(2);
t166 = t142 * qJDD(2);
t151 = qJD(4) * t125 + t145 * t168 + t166;
t171 = qJD(2) * t142;
t162 = qJD(4) * t171;
t134 = qJDD(3) + qJDD(4);
t179 = t134 * pkin(4);
t97 = t151 * t141 + t144 * t162 - t110 - t179;
t187 = t97 - t161;
t140 = sin(qJ(5));
t136 = t140 ^ 2;
t143 = cos(qJ(5));
t186 = (-t143 ^ 2 + t136) * MDP(10);
t135 = qJD(3) + qJD(4);
t133 = t135 ^ 2;
t158 = -g(1) * t111 - g(2) * t112 + t141 * t162;
t118 = t141 * t142 - t144 * t145;
t114 = t118 * qJD(2);
t129 = t141 * pkin(3) + pkin(7);
t130 = -t144 * pkin(3) - pkin(4);
t178 = pkin(3) * qJD(4);
t185 = -qJDD(5) * t129 + (t130 * t135 - t144 * t178 - t114) * qJD(5);
t146 = qJD(5) ^ 2;
t184 = pkin(7) * t146 - t179;
t183 = pkin(3) * t134;
t176 = (t141 * t125 + t144 * t171) * t135;
t119 = t141 * t145 + t144 * t142;
t175 = t119 * qJD(2) * t135;
t173 = qJDD(1) - g(3);
t170 = qJD(5) * t135;
t169 = qJD(5) * t143;
t165 = t135 * t178;
t159 = -t118 * t134 - t133 * t119;
t108 = t144 * t125 - t141 * t171;
t157 = t161 + t176;
t156 = t119 * t146 - t159;
t103 = -t135 * pkin(4) - t108;
t155 = -t134 * pkin(7) - t103 * t135 - t141 * t115 - t151 * t144 + t158;
t153 = -pkin(4) * t170 - pkin(7) * qJDD(5) + qJD(5) * t108;
t99 = t135 * t118;
t152 = qJD(5) * t99 - qJDD(5) * t119 + t118 * t170;
t150 = t129 * t146 + t130 * t134 + t141 * t165 - t175;
t149 = -t151 - t165;
t122 = qJDD(5) * t140 + t146 * t143;
t123 = qJDD(5) * t143 - t146 * t140;
t148 = -0.2e1 * t170 * t186 + t103 * t169 * MDP(15) + t122 * MDP(11) + t123 * MDP(12) + (t136 * MDP(9) + MDP(6)) * t134 + (0.2e1 * t134 * t143 * MDP(10) + t103 * qJD(5) * MDP(14) + 0.2e1 * t135 * t169 * MDP(9) + t187 * MDP(15)) * t140;
t147 = qJD(3) ^ 2;
t117 = t139 * t142 - t177 * t145;
t116 = -t139 * t145 - t177 * t142;
t1 = [-t123 * MDP(14) + t122 * MDP(15) + (MDP(1) + MDP(2)) * t173; (-g(1) * t177 + g(2) * t139 + qJDD(2)) * MDP(2) + (t145 * qJDD(3) - t147 * t142) * MDP(4) + (-qJDD(3) * t142 - t147 * t145) * MDP(5) + t159 * MDP(7) + (-t119 * t134 + t99 * t135) * MDP(8) + (-t156 * MDP(14) + t152 * MDP(15)) * t143 + (t152 * MDP(14) + t156 * MDP(15)) * t140; qJDD(3) * MDP(3) + (g(1) * t117 - g(2) * t116 + t132) * MDP(4) + (-g(1) * t116 - g(2) * t117 - t166) * MDP(5) + (t110 + t161 + t175) * MDP(7) + (-t114 * t135 + t158) * MDP(8) + ((-t162 + t183) * MDP(7) + t149 * MDP(8)) * t144 + (t149 * MDP(7) + (-t115 - t183) * MDP(8)) * t141 + (t185 * MDP(14) + t150 * MDP(15)) * t140 + ((-t150 - t187) * MDP(14) + t185 * MDP(15)) * t143 + t148; (t110 + t157) * MDP(7) + (t108 * t135 + t158) * MDP(8) + (-MDP(7) * t162 - t151 * MDP(8)) * t144 + (-t151 * MDP(7) - t115 * MDP(8)) * t141 + (t153 * MDP(14) + (-t176 + t184) * MDP(15)) * t140 + ((t157 - t97 - t184) * MDP(14) + t153 * MDP(15)) * t143 + t148; qJDD(5) * MDP(13) + t133 * t186 + (t134 * MDP(12) - t173 * MDP(14) + t155 * MDP(15)) * t143 + (-t133 * t143 * MDP(9) + t134 * MDP(11) + t155 * MDP(14) + t173 * MDP(15)) * t140;];
tau = t1;
