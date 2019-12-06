% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:42
% EndTime: 2019-12-05 15:28:44
% DurationCPUTime: 0.57s
% Computational Cost: add. (593->136), mult. (1487->184), div. (0->0), fcn. (1018->4), ass. (0->64)
t143 = sin(pkin(8));
t144 = cos(pkin(8));
t175 = qJD(2) * (t143 ^ 2 + t144 ^ 2);
t167 = qJD(2) * t144;
t126 = qJ(3) * t167 + t143 * qJD(1);
t140 = t144 * qJD(1);
t174 = ((-qJ(3) * qJD(2) * t143 + t140) * t143 - t126 * t144) * MDP(8);
t145 = sin(qJ(4));
t146 = cos(qJ(4));
t124 = t143 * t146 + t144 * t145;
t119 = t124 * qJD(2);
t161 = MDP(14) + MDP(16);
t173 = t119 ^ 2;
t172 = pkin(6) + qJ(3);
t115 = pkin(6) * t167 + t126;
t171 = t115 * t145;
t128 = t172 * t143;
t114 = -qJD(2) * t128 + t140;
t99 = t114 * t146 - t171;
t170 = qJD(5) - t99;
t122 = t124 * qJD(4);
t113 = qJD(2) * t122;
t165 = qJD(2) * t146;
t157 = t144 * t165;
t166 = qJD(2) * t145;
t159 = t143 * t166;
t117 = -t157 + t159;
t164 = qJD(4) * t145;
t158 = t143 * t164;
t163 = qJD(4) * t146;
t121 = -t144 * t163 + t158;
t169 = -t124 * t113 + t121 * t117;
t162 = qJD(2) * qJD(3);
t160 = MDP(15) - MDP(18);
t138 = -pkin(3) * t144 - pkin(2);
t156 = t145 * t162;
t155 = t146 * t162;
t154 = t99 + t171;
t91 = t114 * t164 + t115 * t163 + t143 * t155 + t144 * t156;
t133 = qJD(4) * t157;
t112 = qJD(2) * t158 - t133;
t153 = pkin(4) * t113 + qJ(5) * t112;
t123 = t143 * t145 - t146 * t144;
t152 = -t112 * t123 + t119 * t122;
t100 = t114 * t145 + t115 * t146;
t129 = t172 * t144;
t150 = -t128 * t146 - t129 * t145;
t108 = -t128 * t145 + t129 * t146;
t127 = t138 * qJD(2) + qJD(3);
t95 = pkin(4) * t117 - qJ(5) * t119 + t127;
t149 = -t95 * t119 - t91;
t148 = -t114 * t163 + t143 * t156 - t144 * t155;
t116 = t117 ^ 2;
t104 = pkin(4) * t123 - qJ(5) * t124 + t138;
t103 = pkin(4) * t119 + qJ(5) * t117;
t102 = t133 + (t117 - t159) * qJD(4);
t98 = t124 * qJD(3) + t108 * qJD(4);
t97 = -t123 * qJD(3) + t150 * qJD(4);
t96 = qJD(4) * qJ(5) + t100;
t94 = -qJD(4) * pkin(4) + t170;
t93 = pkin(4) * t122 + qJ(5) * t121 - qJD(5) * t124;
t92 = -qJD(5) * t119 + t153;
t90 = (qJD(5) - t171) * qJD(4) - t148;
t1 = [(t152 + t169) * MDP(17) + (-t121 * t96 + t122 * t94 + t123 * t91 + t124 * t90) * MDP(19) + (t160 * t121 - t161 * t122) * qJD(4); (-t112 * t124 - t119 * t121) * MDP(9) + (-t152 + t169) * MDP(10) + (t113 * t138 + t122 * t127) * MDP(14) + (-t112 * t138 - t121 * t127) * MDP(15) + (t104 * t113 + t117 * t93 + t122 * t95 + t123 * t92) * MDP(16) + (-t108 * t113 + t112 * t150 - t117 * t97 + t119 * t98 - t121 * t94 - t122 * t96 - t123 * t90 + t124 * t91) * MDP(17) + (t104 * t112 - t119 * t93 + t121 * t95 - t92 * t124) * MDP(18) + (t104 * t92 + t108 * t90 - t150 * t91 + t93 * t95 + t94 * t98 + t96 * t97) * MDP(19) + (-MDP(11) * t121 - MDP(12) * t122 - t160 * t97 - t161 * t98) * qJD(4) + (-t174 + (qJ(3) * MDP(8) + (2 * MDP(7))) * t175) * qJD(3); (-t116 - t173) * MDP(17) + (t117 * t96 + (-qJD(5) - t94) * t119 + t153) * MDP(19) - t160 * (-t133 + (t117 + t159) * qJD(4)) + 0.2e1 * t161 * t119 * qJD(4) + (-MDP(7) * t175 + t174) * qJD(2); (-t116 + t173) * MDP(10) + t102 * MDP(11) + (-t119 * t127 - t91) * MDP(14) + t148 * MDP(15) + t149 * MDP(16) + (pkin(4) * t112 - qJ(5) * t113 - (t100 - t96) * t119) * MDP(17) + (t103 * t119 - t148) * MDP(18) + (-pkin(4) * t91 + qJ(5) * t90 - t100 * t94 - t103 * t95 + t170 * t96) * MDP(19) + (t119 * MDP(9) + t127 * MDP(15) - t103 * MDP(16) + (t94 - t170) * MDP(17) - t95 * MDP(18)) * t117 + ((-t143 * t165 - t144 * t166 + t119) * MDP(12) + t154 * MDP(15) + (0.2e1 * qJD(5) - t154) * MDP(18) + t161 * t100) * qJD(4); t119 * t117 * MDP(16) + t102 * MDP(17) + (-qJD(4) ^ 2 - t173) * MDP(18) + (-qJD(4) * t96 - t149) * MDP(19);];
tauc = t1;
