% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:44:15
% EndTime: 2019-12-05 16:44:18
% DurationCPUTime: 0.94s
% Computational Cost: add. (863->149), mult. (2146->210), div. (0->0), fcn. (1390->4), ass. (0->88)
t182 = cos(qJ(3));
t231 = pkin(6) + pkin(7);
t166 = t231 * t182;
t181 = sin(qJ(3));
t213 = qJD(1) * t181;
t150 = qJD(2) * t166 + t213;
t180 = sin(qJ(4));
t144 = t180 * t150;
t207 = qJD(2) * t231;
t149 = t182 * qJD(1) - t181 * t207;
t228 = qJD(3) * pkin(3);
t147 = t149 + t228;
t230 = cos(qJ(4));
t198 = t230 * t147 - t144;
t161 = t180 * t182 + t230 * t181;
t212 = qJD(2) * t161;
t225 = t212 * qJ(5);
t237 = t225 - t198;
t177 = qJD(3) + qJD(4);
t209 = qJD(2) * qJD(3);
t236 = -0.2e1 * t209;
t235 = MDP(5) * t182;
t234 = (t181 ^ 2 - t182 ^ 2) * MDP(6);
t233 = MDP(10) * t181 + MDP(11) * t182;
t232 = t212 ^ 2;
t205 = t230 * t182;
t195 = qJD(2) * t205;
t211 = qJD(2) * t181;
t204 = t180 * t211;
t152 = -t195 + t204;
t175 = -pkin(3) * t182 - pkin(2);
t164 = t175 * qJD(2);
t138 = pkin(4) * t152 + qJD(5) + t164;
t227 = t138 * t212;
t226 = t152 * qJ(5);
t224 = t164 * t212;
t223 = t180 * t181;
t183 = qJD(3) ^ 2;
t222 = t181 * t183;
t221 = t182 * t183;
t122 = pkin(4) * t177 - t237;
t220 = t122 + t237;
t137 = t177 * t161;
t135 = t137 * qJD(2);
t194 = t177 * t223;
t203 = qJD(4) * t230;
t136 = -qJD(3) * t205 - t182 * t203 + t194;
t219 = -t161 * t135 + t136 * t152;
t218 = t230 * t149 - t144;
t217 = t177 * t195;
t210 = qJD(4) * t180;
t208 = t181 * t228;
t206 = qJD(3) * t231;
t146 = t230 * t150;
t202 = t181 * t209;
t172 = pkin(3) * t202;
t201 = pkin(4) * t135 + t172;
t200 = pkin(2) * t236;
t142 = t149 * qJD(3);
t143 = (-t182 * t207 - t213) * qJD(3);
t199 = -t180 * t142 + t230 * t143;
t197 = -t149 * t180 - t146;
t134 = qJD(2) * t194 - t217;
t160 = -t205 + t223;
t193 = -t134 * t160 + t137 * t212;
t192 = -t180 * t147 - t146;
t165 = t231 * t181;
t191 = t180 * t165 - t230 * t166;
t151 = t152 ^ 2;
t190 = t212 * t152 * MDP(12) + (-t151 + t232) * MDP(13) + (t217 + (t152 - t204) * t177) * MDP(14);
t162 = t181 * t206;
t163 = t182 * t206;
t189 = -t230 * t162 - t180 * t163 - t165 * t203 - t166 * t210;
t188 = t192 * qJD(4) + t199;
t187 = t191 * qJD(4) + t180 * t162 - t230 * t163;
t186 = t230 * t142 + t180 * t143 + t147 * t203 - t150 * t210;
t185 = t164 * t152 - t186;
t174 = t230 * pkin(3) + pkin(4);
t133 = -t160 * qJ(5) - t191;
t132 = -t161 * qJ(5) - t230 * t165 - t180 * t166;
t126 = -t225 + t218;
t125 = t197 + t226;
t124 = -t192 - t226;
t121 = t136 * qJ(5) - t161 * qJD(5) + t187;
t120 = -qJ(5) * t137 - qJD(5) * t160 + t189;
t119 = t134 * qJ(5) - qJD(5) * t212 + t188;
t118 = -t135 * qJ(5) - t152 * qJD(5) + t186;
t1 = [(t193 + t219) * MDP(19) + (t118 * t161 - t119 * t160 - t122 * t137 - t124 * t136) * MDP(20) - t233 * t183 + (-t137 * MDP(17) + t136 * MDP(18)) * t177; 0.2e1 * t202 * t235 + t234 * t236 + MDP(7) * t221 - MDP(8) * t222 + (-pkin(6) * t221 + t181 * t200) * MDP(10) + (pkin(6) * t222 + t182 * t200) * MDP(11) + (-t134 * t161 - t136 * t212) * MDP(12) + (-t193 + t219) * MDP(13) + (t175 * t135 + t164 * t137 + t152 * t208 + t172 * t160) * MDP(17) + (-t175 * t134 - t164 * t136 + 0.2e1 * t212 * t208) * MDP(18) + (-t118 * t160 - t119 * t161 - t120 * t152 - t121 * t212 + t122 * t136 - t124 * t137 + t132 * t134 - t133 * t135) * MDP(19) + (t118 * t133 + t124 * t120 + t119 * t132 + t122 * t121 + t201 * (pkin(4) * t160 + t175) + t138 * (pkin(4) * t137 + t208)) * MDP(20) + (-t136 * MDP(14) - t137 * MDP(15) + t187 * MDP(17) - t189 * MDP(18)) * t177; (-pkin(3) * t152 * t211 - t224 - t197 * t177 + (-t146 + (-pkin(3) * t177 - t147) * t180) * qJD(4) + t199) * MDP(17) + (t218 * t177 + (-t177 * t203 - t211 * t212) * pkin(3) + t185) * MDP(18) + (t174 * t134 + (t124 + t125) * t212 + (-t122 + t126) * t152 + (-t135 * t180 + (-t230 * t152 + t180 * t212) * qJD(4)) * pkin(3)) * MDP(19) + (-pkin(4) * t227 + t119 * t174 - t122 * t125 - t124 * t126 + (-t138 * t211 + t118 * t180 + (-t122 * t180 + t230 * t124) * qJD(4)) * pkin(3)) * MDP(20) + t190 + (t233 * pkin(2) - t181 * t235 + t234) * qJD(2) ^ 2; (-t192 * t177 + t188 - t224) * MDP(17) + (t198 * t177 + t185) * MDP(18) + (pkin(4) * t134 - t220 * t152) * MDP(19) + (t220 * t124 + (t119 - t227) * pkin(4)) * MDP(20) + t190; (-t151 - t232) * MDP(19) + (t122 * t212 + t124 * t152 + t201) * MDP(20);];
tauc = t1;
