% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRPRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:54:45
% EndTime: 2019-12-05 15:54:50
% DurationCPUTime: 1.43s
% Computational Cost: add. (869->176), mult. (2303->248), div. (0->0), fcn. (1764->8), ass. (0->90)
t203 = cos(qJ(2));
t234 = qJD(1) * t203;
t219 = qJD(3) - t234;
t197 = cos(pkin(9));
t202 = cos(qJ(4));
t231 = qJD(2) * t202;
t227 = t197 * t231;
t196 = sin(pkin(9));
t199 = sin(qJ(4));
t232 = qJD(2) * t199;
t228 = t196 * t232;
t173 = -t227 + t228;
t201 = cos(qJ(5));
t164 = t201 * t173;
t186 = qJD(4) * t227;
t169 = -t228 * qJD(4) + t186;
t181 = t196 * t202 + t197 * t199;
t177 = t181 * qJD(4);
t170 = qJD(2) * t177;
t175 = qJD(2) * t181;
t198 = sin(qJ(5));
t230 = qJD(5) * t198;
t133 = -qJD(5) * t164 + t201 * t169 - t198 * t170 - t175 * t230;
t148 = t175 * t198 + t164;
t213 = t173 * t198 - t201 * t175;
t206 = qJD(5) * t213 - t169 * t198 - t201 * t170;
t195 = qJD(4) + qJD(5);
t243 = t148 * t195;
t244 = t213 * t195;
t258 = (-t148 ^ 2 + t213 ^ 2) * MDP(17) - t148 * MDP(16) * t213 + (t133 + t243) * MDP(18) + (t206 - t244) * MDP(19);
t182 = (qJD(3) + t234) * qJD(2);
t210 = t181 * t182;
t200 = sin(qJ(2));
t235 = qJD(1) * t200;
t187 = qJD(2) * qJ(3) + t235;
t225 = pkin(6) * qJD(2) + t187;
t167 = t225 * t196;
t168 = t225 * t197;
t214 = t167 * t199 - t168 * t202;
t132 = -pkin(7) * t169 + qJD(4) * t214 - t210;
t142 = -pkin(7) * t173 - t214;
t191 = -pkin(3) * t197 - pkin(2);
t178 = qJD(2) * t191 + t219;
t154 = pkin(4) * t173 + t178;
t257 = t154 * t148 + t142 * t230 + (-t142 * t195 - t132) * t198;
t251 = -t202 * t167 - t168 * t199;
t239 = t202 * t197;
t180 = t196 * t199 - t239;
t253 = t180 * t182;
t131 = -pkin(7) * t170 + t251 * qJD(4) - t253;
t256 = -t198 * t131 + t201 * t132 + t154 * t213;
t237 = t196 ^ 2 + t197 ^ 2;
t254 = MDP(7) * t237;
t252 = qJD(2) * t187 * t237;
t246 = pkin(6) + qJ(3);
t183 = t246 * t196;
t184 = t246 * t197;
t212 = t183 * t199 - t184 * t202;
t250 = qJD(4) * t212 - t181 * t219;
t176 = t180 * qJD(4);
t249 = qJD(5) - t195;
t209 = t180 * t203;
t241 = t183 * t202;
t248 = -qJD(1) * t209 + (qJD(3) * t196 + qJD(4) * t184) * t199 - qJD(3) * t239 + qJD(4) * t241;
t247 = pkin(4) * t175;
t245 = qJD(2) * pkin(2);
t240 = t201 * t142;
t141 = -pkin(7) * t175 + t251;
t138 = qJD(4) * pkin(4) + t141;
t226 = -pkin(4) * t195 - t138;
t224 = t237 * t203;
t223 = t237 * t182;
t220 = t237 * qJD(3);
t218 = pkin(7) * t177 - qJD(5) * (-pkin(7) * t181 - t184 * t199 - t241) + t248;
t217 = -pkin(7) * t176 + qJD(5) * (-pkin(7) * t180 - t212) - t250;
t216 = pkin(4) * t177 - t235;
t152 = t180 * t201 + t181 * t198;
t153 = -t180 * t198 + t181 * t201;
t171 = t181 * t200;
t204 = qJD(2) ^ 2;
t192 = qJD(2) * t235;
t185 = t219 - t245;
t172 = t180 * t200;
t162 = pkin(4) * t180 + t191;
t155 = pkin(4) * t170 + t192;
t146 = -t175 * t203 + t200 * t176;
t145 = -qJD(2) * t209 - qJD(4) * t171;
t136 = qJD(5) * t153 - t176 * t198 + t201 * t177;
t135 = -qJD(5) * t152 - t176 * t201 - t177 * t198;
t1 = [((-t145 * t198 + t146 * t201) * MDP(21) - (t145 * t201 + t146 * t198) * MDP(22) + ((t171 * t198 + t172 * t201) * MDP(21) - (-t171 * t201 + t172 * t198) * MDP(22)) * qJD(5)) * t195 + (MDP(14) * t146 - MDP(15) * t145) * qJD(4) + (-t170 * MDP(14) - t169 * MDP(15) + t206 * MDP(21) - t133 * MDP(22) + MDP(8) * t252 + (-MDP(4) + t254) * t204) * t203 + (MDP(8) * t223 + (-MDP(5) * t197 + MDP(6) * t196 - MDP(3)) * t204 + ((t185 - t234) * MDP(8) + t173 * MDP(14) + t175 * MDP(15) + t148 * MDP(21) - t213 * MDP(22)) * qJD(2)) * t200; (t223 + (-qJD(1) * t224 + t220) * qJD(2)) * MDP(7) + (t187 * t220 + qJ(3) * t223 + ((-t185 - t245) * t200 - t187 * t224) * qJD(1)) * MDP(8) + (t169 * t181 - t175 * t176) * MDP(9) + (-t169 * t180 - t170 * t181 + t173 * t176 - t175 * t177) * MDP(10) + (t191 * t170 + t178 * t177 + (qJD(2) * t180 - t173) * t235) * MDP(14) + (t191 * t169 - t178 * t176) * MDP(15) + (t133 * t153 - t135 * t213) * MDP(16) + (-t133 * t152 - t135 * t148 + t136 * t213 + t153 * t206) * MDP(17) + (t154 * t136 + t216 * t148 + t155 * t152 - t162 * t206) * MDP(21) + (t162 * t133 + t154 * t135 + t155 * t153 - t213 * t216) * MDP(22) + (t135 * MDP(18) - t136 * MDP(19) + (t198 * t218 - t201 * t217) * MDP(21) + (t198 * t217 + t201 * t218) * MDP(22)) * t195 + (-t176 * MDP(11) - t177 * MDP(12) + MDP(14) * t250 + MDP(15) * t248) * qJD(4); (t192 - t252) * MDP(8) + t186 * MDP(15) + (-t206 - t244) * MDP(21) + (t133 - t243) * MDP(22) - t204 * t254 + ((t196 * t231 + t197 * t232 + t175) * MDP(14) + (-t173 - t228) * MDP(15)) * qJD(4); t175 * t173 * MDP(9) + (-t173 ^ 2 + t175 ^ 2) * MDP(10) + (t186 + (t173 - t228) * qJD(4)) * MDP(11) + (-t178 * t175 - t210) * MDP(14) + (t178 * t173 + t253) * MDP(15) + (-(-t141 * t198 - t240) * t195 - t148 * t247 + (t198 * t226 - t240) * qJD(5) + t256) * MDP(21) + (t213 * t247 + (t226 * qJD(5) + t141 * t195 - t131) * t201 + t257) * MDP(22) + t258; (t249 * (-t138 * t198 - t240) + t256) * MDP(21) + ((-t249 * t138 - t131) * t201 + t257) * MDP(22) + t258;];
tauc = t1;
