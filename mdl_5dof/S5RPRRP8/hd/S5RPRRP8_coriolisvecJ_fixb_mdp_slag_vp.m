% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRRP8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:31
% EndTime: 2019-12-31 18:47:33
% DurationCPUTime: 0.79s
% Computational Cost: add. (970->149), mult. (1503->190), div. (0->0), fcn. (636->4), ass. (0->78)
t147 = sin(qJ(4));
t145 = t147 ^ 2;
t149 = cos(qJ(4));
t146 = t149 ^ 2;
t176 = qJD(1) - qJD(3);
t210 = (-MDP(9) + (t145 + t146) * MDP(18)) * t176;
t212 = t176 ^ 2;
t148 = sin(qJ(3));
t184 = qJD(1) * qJ(2);
t151 = -pkin(1) - pkin(2);
t137 = t151 * qJD(1) + qJD(2);
t150 = cos(qJ(3));
t191 = t150 * t137;
t127 = -t148 * t184 + t191;
t136 = -t149 * pkin(4) - t147 * qJ(5) - pkin(3);
t194 = t136 * t176;
t111 = -t127 - t194;
t211 = t176 * t111;
t209 = qJD(4) * t176;
t207 = 0.2e1 * t150 * t209;
t174 = qJD(3) * t184;
t178 = qJD(1) * qJD(2);
t183 = qJD(3) * t148;
t119 = t137 * t183 + t148 * t178 + t150 * t174;
t162 = t150 * qJ(2) + t148 * t151;
t126 = t148 * qJD(2) + t162 * qJD(3);
t206 = t126 * t176 + t119;
t205 = -t148 * qJ(2) + t150 * t151;
t204 = (t145 - t146) * MDP(11);
t177 = MDP(15) + MDP(17);
t163 = pkin(4) * t147 - qJ(5) * t149;
t130 = t163 * qJD(4) - t147 * qJD(5);
t202 = 0.2e1 * qJD(2);
t152 = qJD(4) ^ 2;
t201 = pkin(7) * t152;
t200 = t176 * pkin(3);
t199 = qJD(4) * pkin(4);
t128 = t148 * t137 + t150 * t184;
t197 = t128 * t176;
t196 = t130 * t176;
t135 = -pkin(7) + t162;
t195 = t135 * t152;
t122 = -pkin(7) * t176 + t128;
t193 = t147 * t122;
t192 = t149 * t122;
t188 = qJD(3) * t191 + t150 * t178;
t155 = t148 * t174 - t188;
t115 = t147 * t155;
t107 = qJD(4) * t192 - t115;
t185 = MDP(10) * t149;
t182 = qJD(4) * qJ(5);
t180 = t204 * t209;
t175 = t176 * t185;
t121 = -t127 + t200;
t172 = t121 + t200;
t171 = -t119 - t201;
t170 = pkin(7) * MDP(20) + MDP(18);
t169 = t111 - t194;
t168 = qJD(5) + t193;
t166 = t135 * MDP(20) - MDP(18);
t124 = -t136 - t205;
t125 = t150 * qJD(2) + t205 * qJD(3);
t165 = -t124 * t176 - t111 - t125;
t164 = -(pkin(3) - t205) * t176 - t121 - t125;
t116 = t149 * t155;
t106 = -t116 + (qJD(5) - t193) * qJD(4);
t160 = -t152 * MDP(12) - t106 * MDP(18);
t159 = -t152 * MDP(13) + t107 * MDP(18);
t110 = t119 - t196;
t158 = -t110 + t196 - t201;
t112 = t126 - t130;
t157 = t112 * t176 + t110 - t195;
t156 = t195 - t206;
t154 = t148 * (-t152 - t212);
t131 = t163 * t176;
t114 = t182 + t192;
t113 = t168 - t199;
t1 = [t206 * MDP(8) + t188 * MDP(9) - 0.2e1 * t180 + (t110 * t124 + t111 * t112) * MDP(20) - t125 * t210 + (MDP(5) * t202 + (MDP(6) * t202 - MDP(9) * t183) * qJ(2)) * qJD(1) + (-t156 * MDP(15) + t157 * MDP(17) + (t106 * t135 + t114 * t125) * MDP(20) + (t164 * MDP(16) - t165 * MDP(19) + t166 * t113) * qJD(4) + t160) * t149 + (t156 * MDP(16) + t157 * MDP(19) + (t107 * t135 + t113 * t125) * MDP(20) + (t164 * MDP(15) + t165 * MDP(17) - t166 * t114 + 0.2e1 * t175) * qJD(4) - t159) * t147; (-qJ(2) * MDP(6) - MDP(5)) * qJD(1) ^ 2 + t177 * (t147 * t207 + t149 * t154) + (-MDP(16) + MDP(19)) * (t147 * t154 - t149 * t207) + ((t106 * t149 + t107 * t147 - t211 + (t113 * t149 - t114 * t147) * qJD(4)) * MDP(20) - MDP(8) * t212) * t148 + (-t110 * MDP(20) + ((-t113 * t147 - t114 * t149) * MDP(20) + t210) * t176) * t150; (-t119 - t197) * MDP(8) + t155 * MDP(9) + 0.2e1 * t180 + (t110 * t136 + (-t128 + t130) * t111) * MDP(20) + t127 * t210 + (t171 * MDP(15) + t158 * MDP(17) + (pkin(7) * t106 - t114 * t127) * MDP(20) + ((t127 + t172) * MDP(16) + (-t127 - t169) * MDP(19) + t170 * t113) * qJD(4) - t160) * t149 + ((-t171 + t197) * MDP(16) + (t158 - t197) * MDP(19) + (pkin(7) * t107 - t113 * t127) * MDP(20) + (t172 * MDP(15) + t169 * MDP(17) - t170 * t114 - 0.2e1 * t175) * qJD(4) + t159) * t147 + t177 * (t147 * t127 * qJD(4) - t149 * t197); t116 * MDP(16) + (0.2e1 * qJD(4) * qJD(5) - t116) * MDP(19) + (-t107 * pkin(4) + t106 * qJ(5) + t111 * t131 - t113 * t192 + t168 * t114) * MDP(20) + t177 * t115 + (-t147 * t185 + t204) * t212 - ((-t121 * MDP(15) - t111 * MDP(17) + (t114 - t182) * MDP(18) - t131 * MDP(19)) * t147 + (-t121 * MDP(16) - t131 * MDP(17) + (qJD(5) - t113 - t199) * MDP(18) + t111 * MDP(19)) * t149) * t176; -t147 * t212 * t149 * MDP(17) + (-t145 * t212 - t152) * MDP(19) + (-t114 * qJD(4) - t147 * t211 + t107) * MDP(20);];
tauc = t1;
