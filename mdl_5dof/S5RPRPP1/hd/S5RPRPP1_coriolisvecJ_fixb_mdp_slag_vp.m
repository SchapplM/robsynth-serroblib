% Calculate Coriolis joint torque vector for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:14
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:14:06
% EndTime: 2021-01-15 11:14:10
% DurationCPUTime: 0.96s
% Computational Cost: add. (1049->174), mult. (2481->231), div. (0->0), fcn. (1522->6), ass. (0->79)
t179 = sin(pkin(8));
t181 = cos(pkin(8));
t183 = sin(qJ(3));
t184 = cos(qJ(3));
t163 = t179 * t184 + t181 * t183;
t212 = qJD(1) * t163;
t224 = 0.2e1 * t212;
t173 = sin(pkin(7)) * pkin(1) + pkin(6);
t214 = qJ(4) + t173;
t223 = MDP(14) + MDP(17);
t222 = MDP(6) * (t183 ^ 2 - t184 ^ 2);
t194 = t214 * qJD(1);
t146 = t184 * qJD(2) - t194 * t183;
t204 = t184 * MDP(11);
t220 = -t183 * MDP(10) - t204;
t201 = MDP(12) + MDP(16);
t147 = qJD(2) * t183 + t184 * t194;
t153 = t212 ^ 2;
t202 = qJD(1) * qJD(4);
t138 = t146 * qJD(3) + t184 * t202;
t187 = -t147 * qJD(3) - t183 * t202;
t115 = t138 * t179 - t181 * t187;
t161 = t214 * t184;
t196 = t214 * t183;
t135 = t161 * t179 + t181 * t196;
t219 = t115 * t135;
t215 = t181 * t184;
t162 = t179 * t183 - t215;
t218 = t115 * t162;
t175 = -cos(pkin(7)) * pkin(1) - pkin(2);
t164 = -pkin(3) * t184 + t175;
t211 = qJD(1) * t164;
t154 = qJD(4) + t211;
t208 = qJD(1) * t184;
t198 = t181 * t208;
t209 = qJD(1) * t183;
t155 = t179 * t209 - t198;
t129 = pkin(4) * t155 - qJ(5) * t212 + t154;
t217 = t129 * t212;
t216 = t147 * t179;
t142 = t181 * t147;
t116 = t181 * t138 + t179 * t187;
t144 = qJD(3) * pkin(3) + t146;
t122 = t179 * t144 + t142;
t166 = qJD(1) * t175;
t206 = qJD(3) * t183;
t126 = t146 * t181 - t216;
t203 = qJD(5) - t126;
t200 = MDP(13) - MDP(18);
t199 = pkin(3) * t209;
t197 = qJD(1) * t206;
t195 = qJD(3) * t214;
t121 = t144 * t181 - t216;
t157 = t163 * qJD(3);
t149 = qJD(1) * t157;
t167 = t179 * t197;
t150 = qJD(3) * t198 - t167;
t170 = pkin(3) * t197;
t192 = pkin(4) * t149 - qJ(5) * t150 + t170;
t125 = t146 * t179 + t142;
t191 = qJD(3) * t125 - t115;
t189 = -qJD(4) * t183 - t184 * t195;
t148 = qJD(4) * t184 - t183 * t195;
t127 = t148 * t179 - t181 * t189;
t128 = t181 * t148 + t179 * t189;
t136 = t181 * t161 - t179 * t196;
t188 = t115 * t163 + t127 * t212 - t128 * t155 + t135 * t150 - t136 * t149;
t185 = qJD(3) ^ 2;
t174 = -pkin(3) * t181 - pkin(4);
t171 = pkin(3) * t179 + qJ(5);
t160 = qJD(3) * t215 - t179 * t206;
t131 = pkin(4) * t162 - qJ(5) * t163 + t164;
t130 = pkin(4) * t212 + qJ(5) * t155 + t199;
t123 = pkin(3) * t206 + pkin(4) * t157 - qJ(5) * t160 - qJD(5) * t163;
t120 = qJD(3) * qJ(5) + t122;
t119 = -qJD(3) * pkin(4) + qJD(5) - t121;
t118 = -qJD(5) * t212 + t192;
t114 = qJD(3) * qJD(5) + t116;
t1 = [(t149 * t164 + t154 * t157) * MDP(12) + (t150 * t164 + t154 * t160) * MDP(13) + (-t116 * t162 - t121 * t160 - t122 * t157 + t188) * MDP(14) + (t116 * t136 - t121 * t127 + t122 * t128 + t219) * MDP(15) + (t118 * t162 + t123 * t155 + t129 * t157 + t131 * t149) * MDP(16) + (-t114 * t162 + t119 * t160 - t120 * t157 + t188) * MDP(17) + (-t118 * t163 - t123 * t212 - t129 * t160 - t131 * t150) * MDP(18) + (t114 * t136 + t118 * t131 + t119 * t127 + t120 * t128 + t123 * t129 + t219) * MDP(19) + (t184 * MDP(7) - t183 * MDP(8) + (-MDP(10) * t184 + MDP(11) * t183) * t173) * t185 + (t166 * t204 - t200 * t128 - t201 * t127 + (t175 * t204 - 0.2e1 * t222) * qJD(1) + (0.2e1 * MDP(5) * t208 + 0.2e1 * t166 * MDP(10) + ((qJD(1) * t162 + t155) * MDP(12) + MDP(13) * t224 + (t154 + t211) * MDP(15)) * pkin(3)) * t183) * qJD(3); (t116 * t163 - t121 * t157 + t122 * t160 + t218) * MDP(15) + (t114 * t163 + t119 * t157 + t120 * t160 + t218) * MDP(19) + t220 * t185 + (-t157 * t201 - t160 * t200) * qJD(3) + t223 * (-t163 * t149 + t150 * t162 - t160 * t155 + t157 * t212); (-t154 * t212 - t155 * t199 + t191) * MDP(12) + (qJD(3) * t126 + t154 * t155 - t199 * t212 - t116) * MDP(13) + ((t122 - t125) * t212 + (-t121 + t126) * t155 + (-t149 * t179 - t150 * t181) * pkin(3)) * MDP(14) + (t121 * t125 - t122 * t126 + (-t115 * t181 + t116 * t179 - t154 * t209) * pkin(3)) * MDP(15) + (-t130 * t155 + t191 - t217) * MDP(16) + (-t149 * t171 + t150 * t174 + (t120 - t125) * t212 + (t119 - t203) * t155) * MDP(17) + (-t129 * t155 + t130 * t212 + (0.2e1 * qJD(5) - t126) * qJD(3) + t116) * MDP(18) + (t114 * t171 + t115 * t174 - t119 * t125 + t120 * t203 - t129 * t130) * MDP(19) + (t220 * t166 + (-t183 * t184 * MDP(5) + t222) * qJD(1)) * qJD(1); (t121 * t212 + t122 * t155 + t170) * MDP(15) + (t120 * t155 + (-qJD(5) - t119) * t212 + t192) * MDP(19) + t200 * (-t167 + (-t155 + t198) * qJD(3)) + t223 * (-t155 ^ 2 - t153) + t201 * qJD(3) * t224; t212 * t155 * MDP(16) + (-t167 + (t155 + t198) * qJD(3)) * MDP(17) + (-t153 - t185) * MDP(18) + (-qJD(3) * t120 + t115 + t217) * MDP(19);];
tauc = t1;
