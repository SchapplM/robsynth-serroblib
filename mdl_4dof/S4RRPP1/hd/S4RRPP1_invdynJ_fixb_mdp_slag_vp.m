% Calculate vector of inverse dynamics joint torques for
% S4RRPP1
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
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RRPP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:33:10
% EndTime: 2019-03-08 18:33:11
% DurationCPUTime: 0.34s
% Computational Cost: add. (337->107), mult. (509->133), div. (0->0), fcn. (265->10), ass. (0->66)
t146 = cos(qJ(2));
t170 = pkin(1) * qJD(2);
t156 = qJD(1) * t170;
t144 = sin(qJ(2));
t160 = qJDD(1) * t144;
t176 = pkin(1) * t160 + t146 * t156;
t141 = qJ(1) + qJ(2);
t135 = sin(t141);
t125 = g(1) * t135;
t136 = cos(t141);
t175 = -g(2) * t136 + t125;
t174 = pkin(1) * t146;
t145 = sin(qJ(1));
t173 = g(1) * t145;
t171 = pkin(1) * qJD(1);
t142 = sin(pkin(6));
t169 = MDP(7) * t142;
t128 = pkin(2) + t174;
t143 = cos(pkin(6));
t168 = t128 * t143;
t167 = t142 * t144;
t166 = t143 * t144;
t165 = t143 * t146;
t139 = qJD(1) + qJD(2);
t157 = t146 * t171;
t107 = pkin(2) * t139 + t157;
t158 = t144 * t171;
t113 = t143 * t158;
t96 = t142 * t107 + t113;
t164 = pkin(1) * t166 + t142 * t128;
t163 = MDP(7) + MDP(10);
t162 = -qJD(1) - t139;
t159 = pkin(1) * t167;
t129 = qJDD(1) * t174;
t138 = qJDD(1) + qJDD(2);
t100 = pkin(2) * t138 - t144 * t156 + t129;
t91 = t143 * t100 - t176 * t142;
t92 = t142 * t100 + t176 * t143;
t134 = pkin(6) + t141;
t121 = sin(t134);
t122 = cos(t134);
t127 = pkin(2) * t136;
t154 = t122 * pkin(3) + t121 * qJ(4) + t127;
t153 = -pkin(3) * t121 + t122 * qJ(4);
t152 = t138 * qJ(4) + t92;
t90 = -t138 * pkin(3) + qJDD(4) - t91;
t95 = t107 * t143 - t142 * t158;
t150 = -g(1) * t122 - g(2) * t121 + t152;
t149 = -g(1) * t121 + g(2) * t122 + t90;
t148 = t138 * MDP(4) + (g(1) * t136 + g(2) * t135) * MDP(6) + (t129 + t175) * MDP(5);
t147 = cos(qJ(1));
t137 = t147 * pkin(1);
t132 = t139 * qJD(4);
t123 = -pkin(2) * t143 - pkin(3);
t120 = pkin(2) * t142 + qJ(4);
t114 = t165 * t170;
t105 = (t165 - t167) * t171;
t104 = (t142 * t146 + t166) * t170;
t103 = t142 * t157 + t113;
t102 = -pkin(3) + t159 - t168;
t101 = qJ(4) + t164;
t99 = -qJD(2) * t159 + qJD(4) + t114;
t94 = qJ(4) * t139 + t96;
t93 = -pkin(3) * t139 + qJD(4) - t95;
t89 = t132 + t152;
t1 = [qJDD(1) * MDP(1) + (-g(2) * t147 + t173) * MDP(2) + (g(1) * t147 + g(2) * t145) * MDP(3) + (t92 * t164 + t96 * t114 + t91 * t168 - t95 * t104 + pkin(2) * t125 - g(2) * (t127 + t137)) * MDP(7) + (-t102 * t138 - t104 * t139 - t149) * MDP(8) + (t101 * t138 + t139 * t99 + t132 + t150) * MDP(9) + (t89 * t101 + t94 * t99 + t90 * t102 + t93 * t104 - g(1) * (-pkin(2) * t135 + t153) - g(2) * (t137 + t154)) * MDP(10) + (t163 * t173 + (t162 * MDP(6) * qJD(2) + t138 * MDP(5)) * t146 + ((-qJDD(1) - t138) * MDP(6) - t91 * t169 + (t162 * MDP(5) - t96 * t169) * qJD(2)) * t144) * pkin(1) + t148; (t103 * t95 - t105 * t96) * MDP(7) + (t103 * t139 - t123 * t138 - t149) * MDP(8) + (-t105 * t139 + t120 * t138 + 0.2e1 * t132 + t150) * MDP(9) + (-g(1) * t153 - g(2) * t154 - t93 * t103 + t89 * t120 + t90 * t123 + (qJD(4) - t105) * t94) * MDP(10) + ((t142 * t92 + t143 * t91 + t175) * MDP(7) + MDP(10) * t125) * pkin(2) + (-MDP(6) * t160 + (t144 * MDP(5) + t146 * MDP(6)) * qJD(1) * (-qJD(2) + t139)) * pkin(1) + t148; t163 * (qJDD(3) - g(3)); -t138 * MDP(8) - t139 ^ 2 * MDP(9) + (-t139 * t94 + t149) * MDP(10);];
tau  = t1;
