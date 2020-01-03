% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:56
% EndTime: 2019-12-31 17:04:58
% DurationCPUTime: 0.96s
% Computational Cost: add. (740->153), mult. (2023->229), div. (0->0), fcn. (1431->6), ass. (0->94)
t216 = sin(pkin(7));
t217 = cos(pkin(7));
t219 = sin(qJ(2));
t221 = cos(qJ(2));
t200 = -t216 * t219 + t217 * t221;
t190 = t200 * qJD(1);
t220 = cos(qJ(4));
t181 = t220 * t190;
t201 = t216 * t221 + t217 * t219;
t191 = t201 * qJD(2);
t184 = qJD(1) * t191;
t237 = qJD(1) * qJD(2);
t233 = t221 * t237;
t234 = t219 * t237;
t185 = -t216 * t234 + t217 * t233;
t192 = t201 * qJD(1);
t218 = sin(qJ(4));
t240 = qJD(4) * t218;
t131 = qJD(4) * t181 - t218 * t184 + t220 * t185 - t192 * t240;
t226 = t190 * t218 + t220 * t192;
t132 = qJD(4) * t226 + t220 * t184 + t185 * t218;
t154 = -t192 * t218 + t181;
t213 = qJD(2) + qJD(4);
t246 = t154 * t213;
t247 = t226 * t213;
t261 = (-t132 + t247) * MDP(16) + (-t154 ^ 2 + t226 ^ 2) * MDP(14) - t154 * t226 * MDP(13) + (t131 - t246) * MDP(15);
t249 = -qJ(3) - pkin(5);
t208 = t249 * t219;
t205 = qJD(1) * t208;
t248 = qJD(2) * pkin(2);
t199 = t205 + t248;
t209 = t249 * t221;
t206 = qJD(1) * t209;
t245 = t217 * t206;
t157 = t216 * t199 - t245;
t251 = pkin(6) * t190;
t142 = t157 + t251;
t235 = -pkin(2) * t221 - pkin(1);
t228 = t235 * qJD(1);
t207 = qJD(3) + t228;
t162 = -pkin(3) * t190 + t207;
t260 = t142 * t240 - t162 * t154;
t257 = -0.2e1 * t237;
t256 = MDP(4) * t219;
t255 = MDP(5) * (t219 ^ 2 - t221 ^ 2);
t254 = qJD(4) - t213;
t232 = qJD(2) * t249;
t188 = qJD(3) * t221 + t219 * t232;
t174 = t188 * qJD(1);
t189 = -qJD(3) * t219 + t221 * t232;
t175 = t189 * qJD(1);
t143 = -t174 * t216 + t217 * t175;
t136 = -pkin(6) * t185 + t143;
t144 = t217 * t174 + t216 * t175;
t137 = -pkin(6) * t184 + t144;
t253 = t220 * t136 - t218 * t137 - t162 * t226;
t252 = pkin(2) * t216;
t250 = pkin(6) * t192;
t195 = t216 * t206;
t222 = qJD(2) ^ 2;
t244 = t219 * t222;
t243 = t221 * t222;
t223 = qJD(1) ^ 2;
t242 = t221 * t223;
t148 = t217 * t188 + t216 * t189;
t161 = t217 * t205 + t195;
t165 = t216 * t208 - t217 * t209;
t238 = t219 * qJD(1);
t236 = t219 * t248;
t231 = pkin(1) * t257;
t147 = -t188 * t216 + t217 * t189;
t156 = t217 * t199 + t195;
t160 = -t205 * t216 + t245;
t164 = t217 * t208 + t209 * t216;
t141 = qJD(2) * pkin(3) + t156 - t250;
t227 = -t218 * t141 - t220 * t142;
t225 = t200 * t220 - t201 * t218;
t159 = t200 * t218 + t201 * t220;
t212 = pkin(2) * t217 + pkin(3);
t211 = pkin(2) * t234;
t194 = t200 * qJD(2);
t177 = -pkin(3) * t200 + t235;
t170 = pkin(3) * t191 + t236;
t169 = pkin(2) * t238 + pkin(3) * t192;
t163 = pkin(3) * t184 + t211;
t150 = pkin(6) * t200 + t165;
t149 = -pkin(6) * t201 + t164;
t146 = t161 - t250;
t145 = t160 - t251;
t139 = -pkin(6) * t191 + t148;
t138 = -pkin(6) * t194 + t147;
t134 = qJD(4) * t159 + t220 * t191 + t194 * t218;
t133 = qJD(4) * t225 - t191 * t218 + t194 * t220;
t1 = [0.2e1 * t233 * t256 + t255 * t257 + MDP(6) * t243 - MDP(7) * t244 + (-pkin(5) * t243 + t219 * t231) * MDP(9) + (pkin(5) * t244 + t221 * t231) * MDP(10) + (-t143 * t201 + t144 * t200 - t147 * t192 + t148 * t190 - t156 * t194 - t157 * t191 - t164 * t185 - t165 * t184) * MDP(11) + (t143 * t164 + t144 * t165 + t156 * t147 + t157 * t148 + (t207 + t228) * t236) * MDP(12) + (t131 * t159 + t133 * t226) * MDP(13) + (t131 * t225 - t132 * t159 + t133 * t154 - t134 * t226) * MDP(14) + (t177 * t132 + t162 * t134 - t154 * t170 - t163 * t225) * MDP(18) + (t177 * t131 + t162 * t133 + t163 * t159 + t170 * t226) * MDP(19) + (t133 * MDP(15) - t134 * MDP(16) + (t138 * t220 - t139 * t218 + (-t149 * t218 - t150 * t220) * qJD(4)) * MDP(18) + (-t138 * t218 - t139 * t220 - (t149 * t220 - t150 * t218) * qJD(4)) * MDP(19)) * t213; -t242 * t256 + t223 * t255 + ((t157 + t160) * t192 + (t156 - t161) * t190 + (-t184 * t216 - t185 * t217) * pkin(2)) * MDP(11) + (-t156 * t160 - t157 * t161 + (t143 * t217 + t144 * t216 - t207 * t238) * pkin(2)) * MDP(12) + (t169 * t154 - (t145 * t220 - t146 * t218) * t213 + ((-t212 * t218 - t220 * t252) * t213 + t227) * qJD(4) + t253) * MDP(18) + (-t220 * t137 - t218 * t136 - t169 * t226 + (t145 * t218 + t146 * t220) * t213 + (-(t212 * t220 - t218 * t252) * t213 - t220 * t141) * qJD(4) + t260) * MDP(19) + (MDP(9) * t219 * t223 + MDP(10) * t242) * pkin(1) + t261; (-t190 ^ 2 - t192 ^ 2) * MDP(11) + (t156 * t192 - t157 * t190 + t211) * MDP(12) + (t132 + t247) * MDP(18) + (t131 + t246) * MDP(19); (t254 * t227 + t253) * MDP(18) + ((-t142 * t213 - t136) * t218 + (-t254 * t141 - t137) * t220 + t260) * MDP(19) + t261;];
tauc = t1;
