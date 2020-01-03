% Calculate vector of inverse dynamics joint torques for
% S4RPRR5
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RPRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:42
% EndTime: 2019-12-31 16:51:43
% DurationCPUTime: 0.48s
% Computational Cost: add. (397->103), mult. (598->134), div. (0->0), fcn. (322->6), ass. (0->54)
t131 = cos(qJ(1));
t173 = sin(qJ(1));
t182 = g(1) * t173 - g(2) * t131;
t155 = qJD(1) - qJD(3);
t181 = t155 ^ 2;
t132 = -pkin(1) - pkin(2);
t160 = qJ(2) * qJD(1);
t180 = -qJD(3) * t160 + t132 * qJDD(1) + qJDD(2);
t179 = -qJDD(2) + t182;
t158 = qJD(4) * t155;
t110 = t132 * qJD(1) + qJD(2);
t156 = (qJD(1) * qJD(2));
t157 = (qJ(2) * qJDD(1));
t178 = qJD(3) * t110 + t156 + t157;
t128 = sin(qJ(3));
t130 = cos(qJ(3));
t164 = t130 * qJ(2) + t128 * t132;
t167 = pkin(1) * qJDD(1);
t176 = t167 + t179;
t102 = -t173 * t128 - t131 * t130;
t103 = t128 * t131 - t173 * t130;
t135 = g(1) * t102 + g(2) * t103 + t180 * t128 + t178 * t130;
t139 = -g(1) * t103 + g(2) * t102 + t178 * t128 - t180 * t130;
t171 = pkin(3) * t155;
t98 = t110 * t130 - t128 * t160;
t94 = -t98 + t171;
t175 = (t94 + t98 + t171) * qJD(4) - pkin(6) * qJDD(4);
t127 = sin(qJ(4));
t129 = cos(qJ(4));
t174 = MDP(15) * t129 - MDP(16) * t127 + MDP(8);
t123 = qJDD(1) - qJDD(3);
t172 = pkin(3) * t123;
t169 = t155 * (t110 * t128 + t130 * t160);
t168 = (qJD(2) * t128 + t164 * qJD(3)) * t155;
t165 = t127 * t129;
t125 = t127 ^ 2;
t163 = -t129 ^ 2 + t125;
t149 = -qJ(2) * t128 + t130 * t132;
t133 = qJD(4) ^ 2;
t108 = qJDD(4) * t129 - t127 * t133;
t107 = qJDD(4) * t127 + t129 * t133;
t148 = -t172 - t139;
t146 = g(1) * t131 + g(2) * t173;
t143 = pkin(6) * t123 + t155 * t94 - t135;
t142 = -t146 + (2 * t156);
t104 = pkin(3) - t149;
t105 = -pkin(6) + t164;
t96 = qJD(2) * t130 + t149 * qJD(3);
t140 = -qJDD(4) * t105 + (-t104 * t155 - t94 - t96) * qJD(4);
t138 = pkin(6) * t133 - t148 + t169 + t172;
t137 = -t104 * t123 + t105 * t133 + t148 - t168;
t136 = (-t123 * t125 - 0.2e1 * t158 * t165) * MDP(10) + 0.2e1 * (-t123 * t165 + t163 * t158) * MDP(11) + t107 * MDP(12) + t108 * MDP(13) - t123 * MDP(7);
t134 = qJD(1) ^ 2;
t1 = [qJDD(1) * MDP(1) + t182 * MDP(2) + t146 * MDP(3) + (0.2e1 * t167 + t179) * MDP(4) + (t142 + (2 * t157)) * MDP(5) + (t176 * pkin(1) + (t142 + t157) * qJ(2)) * MDP(6) + (-t149 * t123 + t139 + t168) * MDP(8) + (t164 * t123 + t155 * t96 + t135) * MDP(9) + (t140 * t127 - t137 * t129) * MDP(15) + (t137 * t127 + t140 * t129) * MDP(16) - t136; -qJDD(1) * MDP(4) - t134 * MDP(5) + (-qJ(2) * t134 - t176) * MDP(6) + (-MDP(9) * t181 - t174 * t123 + 0.2e1 * (MDP(15) * t127 + MDP(16) * t129) * t158) * t130 + (-t107 * MDP(15) - t108 * MDP(16) + t123 * MDP(9) - t174 * t181) * t128; (-t139 - t169) * MDP(8) + (-t155 * t98 - t135) * MDP(9) + (-t138 * MDP(15) + t175 * MDP(16)) * t129 + (t175 * MDP(15) + t138 * MDP(16)) * t127 + t136; qJDD(4) * MDP(14) + t163 * MDP(11) * t181 + (-t123 * MDP(13) + g(3) * MDP(15) + t143 * MDP(16)) * t129 + (-MDP(10) * t129 * t181 - t123 * MDP(12) + t143 * MDP(15) - g(3) * MDP(16)) * t127;];
tau = t1;
