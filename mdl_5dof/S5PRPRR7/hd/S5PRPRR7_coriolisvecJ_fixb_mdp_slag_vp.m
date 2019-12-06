% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5PRPRR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:00:48
% EndTime: 2019-12-05 16:00:52
% DurationCPUTime: 0.90s
% Computational Cost: add. (424->137), mult. (957->201), div. (0->0), fcn. (595->6), ass. (0->76)
t141 = sin(qJ(5));
t142 = sin(qJ(4));
t144 = cos(qJ(5));
t145 = cos(qJ(4));
t122 = t141 * t145 + t142 * t144;
t138 = qJD(4) + qJD(5);
t154 = t122 * t138;
t103 = qJD(2) * t154;
t190 = qJD(5) - t138;
t193 = MDP(9) * (t142 ^ 2 - t145 ^ 2);
t173 = qJD(4) * t145;
t192 = qJD(5) * t145 + t173;
t191 = t145 * MDP(14);
t147 = -pkin(2) - pkin(6);
t189 = pkin(7) - t147;
t188 = qJD(2) * pkin(2);
t177 = qJD(2) * qJ(3);
t143 = sin(qJ(2));
t179 = qJD(1) * t143;
t128 = t177 + t179;
t146 = cos(qJ(2));
t187 = t128 * t146;
t178 = qJD(1) * t146;
t160 = qJD(3) - t178;
t121 = qJD(2) * t147 + t160;
t176 = qJD(2) * t142;
t111 = -pkin(7) * t176 + t121 * t142;
t185 = t144 * t111;
t184 = t144 * t145;
t149 = qJD(2) ^ 2;
t183 = t149 * MDP(6);
t170 = qJD(2) * qJD(4);
t165 = t142 * t170;
t172 = qJD(5) * t142;
t167 = t141 * t172;
t182 = -qJD(2) * t167 - t141 * t165;
t148 = qJD(4) ^ 2;
t180 = t148 + t149;
t175 = qJD(2) * t145;
t174 = qJD(4) * t142;
t169 = pkin(4) * t175;
t168 = t143 * t175;
t166 = t144 * t175;
t125 = t189 * t145;
t136 = pkin(4) * t142 + qJ(3);
t112 = -pkin(7) * t175 + t145 * t121;
t110 = qJD(4) * pkin(4) + t112;
t164 = -pkin(4) * t138 - t110;
t163 = t180 * t146;
t162 = -t128 + t179;
t161 = -0.2e1 * t165;
t132 = pkin(4) * t173 + qJD(3);
t159 = t138 * t184;
t123 = -t141 * t142 + t184;
t130 = qJD(1) * t168;
t105 = t130 + (pkin(7) * qJD(2) - t121) * t174;
t106 = t121 * t173 + (-pkin(7) * t173 + t142 * t179) * qJD(2);
t115 = t141 * t176 - t166;
t118 = qJD(2) * t136 + t179;
t157 = t144 * t105 - t141 * t106 + t118 * t115;
t156 = -t162 + t177;
t116 = t122 * qJD(2);
t155 = -t115 * t116 * MDP(15) + (-t182 + (-t115 - t166) * t138) * MDP(18) + (t115 ^ 2 - t116 ^ 2) * MDP(16) + (t116 * t138 - t103) * MDP(17);
t153 = t123 * t138;
t152 = -t141 * t174 - t167;
t151 = t118 * t116 + (t190 * t111 - t105) * t141;
t126 = (qJD(3) + t178) * qJD(2);
t150 = qJD(2) * t160 - t147 * t148 + t126;
t127 = t160 - t188;
t124 = t189 * t142;
t120 = qJD(4) * t125;
t119 = t189 * t174;
t113 = (t132 + t178) * qJD(2);
t108 = t159 + t152;
t104 = qJD(2) * t159 + t182;
t1 = [t187 * qJD(2) * MDP(7) + (0.2e1 * qJD(4) * t168 + t142 * t163) * MDP(13) + t163 * t191 + (-t149 * MDP(4) + t183 + (((t172 + t174) * t144 + t192 * t141) * t138 + qJD(2) * t116) * MDP(20) + (-(-t192 * t144 - t152) * t138 - qJD(2) * t115) * MDP(21)) * t146 + ((-MDP(3) + MDP(5)) * t149 + (t126 + (t127 - t178) * qJD(2)) * MDP(7) + t161 * MDP(14) + (qJD(2) * t153 + t104) * MDP(20) - 0.2e1 * t103 * MDP(21)) * t143; 0.2e1 * qJD(2) * qJD(3) * MDP(6) + (qJ(3) * t126 + qJD(3) * t128 + (-t187 + (-t127 - t188) * t143) * qJD(1)) * MDP(7) + 0.2e1 * t170 * t193 - t148 * t142 * MDP(10) + (t142 * t150 + t156 * t173) * MDP(13) - t156 * t174 * MDP(14) + (-t103 * t123 + t115 * t154) * MDP(15) + (t103 * t122 - t104 * t123 + t108 * t115 + t116 * t154) * MDP(16) + (t132 * t116 + t136 * t104 + t113 * t122 + t118 * t108 + (-t146 * t116 - t143 * t153) * qJD(1)) * MDP(20) + (-t132 * t115 - t136 * t103 + t113 * t123 - t118 * t154 + (t146 * t115 + t143 * t154) * qJD(1)) * MDP(21) + (-t148 * MDP(11) + t150 * MDP(14) + MDP(8) * t161) * t145 + (-t154 * MDP(17) - t108 * MDP(18) + (t119 * t144 + t120 * t141 + (t124 * t144 + t125 * t141) * qJD(5)) * MDP(20) + (-t119 * t141 + t120 * t144 - (t124 * t141 - t125 * t144) * qJD(5)) * MDP(21)) * t138; -t183 + (-MDP(20) * t154 - MDP(21) * t108) * t138 + (-t116 * MDP(20) + t115 * MDP(21) + MDP(7) * t162) * qJD(2) + (-t142 * MDP(13) - t191) * t180; (-t128 * t175 + t130) * MDP(13) - t162 * MDP(14) * t176 + (-(-t112 * t141 - t185) * t138 - t116 * t169 + (t141 * t164 - t185) * qJD(5) + t157) * MDP(20) + (t115 * t169 + (qJD(5) * t164 + t112 * t138 - t106) * t144 + t151) * MDP(21) + t155 + (t145 * t142 * MDP(8) - t193) * t149; (t157 + t190 * (-t110 * t141 - t185)) * MDP(20) + ((-t190 * t110 - t106) * t144 + t151) * MDP(21) + t155;];
tauc = t1;
