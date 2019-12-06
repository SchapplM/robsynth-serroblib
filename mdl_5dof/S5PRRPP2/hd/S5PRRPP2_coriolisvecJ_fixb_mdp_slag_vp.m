% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:10:11
% EndTime: 2019-12-05 16:10:15
% DurationCPUTime: 1.03s
% Computational Cost: add. (950->181), mult. (2395->252), div. (0->0), fcn. (1584->6), ass. (0->91)
t204 = sin(qJ(3));
t206 = cos(qJ(3));
t254 = -qJ(4) - pkin(6);
t221 = qJD(3) * t254;
t172 = qJD(4) * t206 + t204 * t221;
t202 = sin(pkin(8));
t203 = cos(pkin(8));
t213 = -qJD(4) * t204 + t206 * t221;
t183 = t202 * t204 - t203 * t206;
t207 = cos(qJ(2));
t216 = t183 * t207;
t246 = qJD(1) * t216 + t203 * t172 + t202 * t213;
t261 = MDP(12) + MDP(15);
t260 = MDP(5) * t204;
t259 = (t204 ^ 2 - t206 ^ 2) * MDP(6);
t237 = qJD(3) * t206;
t238 = qJD(3) * t204;
t257 = t202 * t238 - t203 * t237;
t239 = qJD(2) * t206;
t240 = qJD(2) * t204;
t256 = -t202 * t239 - t203 * t240;
t205 = sin(qJ(2));
t242 = qJD(1) * t205;
t255 = pkin(3) * t238 - t242;
t184 = t202 * t206 + t203 * t204;
t177 = t184 * qJD(3);
t178 = t184 * qJD(2);
t174 = t178 ^ 2;
t253 = qJD(2) * pkin(2);
t191 = qJD(2) * pkin(6) + t242;
t241 = qJD(1) * t207;
t217 = qJD(4) + t241;
t155 = -t191 * t238 + (-qJ(4) * t238 + t217 * t206) * qJD(2);
t210 = -t191 * t237 + (-qJ(4) * t237 - t217 * t204) * qJD(2);
t132 = t155 * t202 - t203 * t210;
t189 = t254 * t206;
t224 = t254 * t204;
t157 = -t189 * t202 - t203 * t224;
t252 = t132 * t157;
t168 = t184 * t205;
t251 = t132 * t168;
t220 = qJ(4) * qJD(2) + t191;
t171 = t220 * t206;
t250 = t171 * t202;
t160 = t203 * t171;
t208 = qJD(3) ^ 2;
t249 = t204 * t208;
t248 = t206 * t208;
t247 = t172 * t202 - t184 * t241 - t203 * t213;
t133 = t203 * t155 + t202 * t210;
t170 = t220 * t204;
t162 = qJD(3) * pkin(3) - t170;
t142 = t202 * t162 + t160;
t233 = qJD(2) * qJD(3);
t223 = t204 * t233;
t182 = pkin(3) * t223 + qJD(2) * t242;
t230 = -pkin(3) * t206 - pkin(2);
t181 = t230 * qJD(2) + qJD(4) - t241;
t243 = MDP(13) * t181;
t225 = t203 * t239;
t175 = t202 * t240 - t225;
t236 = t175 * MDP(14);
t235 = t206 * MDP(11);
t145 = -t170 * t203 - t250;
t234 = qJD(5) - t145;
t232 = qJD(3) * qJD(5);
t222 = t206 * t233;
t219 = -pkin(4) * t177 - qJ(5) * t257 + qJD(5) * t184 - t255;
t141 = t162 * t203 - t250;
t166 = qJD(2) * t177;
t190 = t202 * t223;
t167 = t203 * t222 - t190;
t215 = pkin(4) * t166 - qJ(5) * t167 + t182;
t212 = -0.2e1 * qJD(3) * t253;
t158 = -t203 * t189 + t202 * t224;
t211 = t132 * t184 + t157 * t167 - t158 * t166 - t246 * t175 + t247 * t178;
t209 = qJD(2) ^ 2;
t198 = -pkin(3) * t203 - pkin(4);
t196 = pkin(3) * t202 + qJ(5);
t169 = t183 * t205;
t156 = pkin(4) * t183 - qJ(5) * t184 + t230;
t148 = pkin(3) * t240 + pkin(4) * t178 + qJ(5) * t175;
t147 = -qJD(2) * t216 - t205 * t177;
t146 = t257 * t205 + t256 * t207;
t144 = -t170 * t202 + t160;
t140 = qJD(3) * qJ(5) + t142;
t139 = pkin(4) * t175 - qJ(5) * t178 + t181;
t138 = -qJD(3) * pkin(4) + qJD(5) - t141;
t134 = -qJD(5) * t178 + t215;
t131 = t232 + t133;
t1 = [(-t133 * t169 + t141 * t146 + t142 * t147 + t251) * MDP(13) + (-t131 * t169 - t138 * t146 + t140 * t147 + t251) * MDP(17) + (MDP(14) * t146 + MDP(16) * t147) * qJD(3) + (-t182 * MDP(13) - t166 * MDP(14) + t167 * MDP(16) - t134 * MDP(17) - t209 * MDP(4) + 0.2e1 * (-t204 * MDP(10) - t235) * t233) * t207 + (-t209 * MDP(3) + (-MDP(16) * t178 + t139 * MDP(17) + t236 + t243) * qJD(2) + (-MDP(10) * t206 + MDP(11) * t204) * (t208 + t209)) * t205 + t261 * (-t146 * t178 - t147 * t175 + t169 * t166 + t167 * t168); 0.2e1 * t222 * t260 - 0.2e1 * t233 * t259 + MDP(7) * t248 - MDP(8) * t249 + (-pkin(6) * t248 + t204 * t212) * MDP(10) + (pkin(6) * t249 + t206 * t212) * MDP(11) + (-t133 * t183 + t141 * t257 - t142 * t177 + t211) * MDP(12) + (t133 * t158 - t247 * t141 + t246 * t142 + t255 * t181 + t182 * t230 + t252) * MDP(13) + (-t247 * qJD(3) + t134 * t183 + t139 * t177 + t156 * t166 - t219 * t175) * MDP(14) + (-t131 * t183 - t138 * t257 - t140 * t177 + t211) * MDP(15) + (t246 * qJD(3) - t134 * t184 + t139 * t257 - t156 * t167 + t219 * t178) * MDP(16) + (t131 * t158 + t134 * t156 + t247 * t138 - t219 * t139 + t246 * t140 + t252) * MDP(17); (t141 * t144 - t142 * t145) * MDP(13) + (qJD(3) * t144 - t132) * MDP(14) + (-t166 * t196 + t167 * t198) * MDP(15) + (-qJD(3) * t145 + t133 + 0.2e1 * t232) * MDP(16) + (t131 * t196 + t132 * t198 - t138 * t144 - t139 * t148 + t234 * t140) * MDP(17) + ((t142 - t144) * MDP(12) - t139 * MDP(14) + (t140 - t144) * MDP(15) + t148 * MDP(16)) * t178 + ((-t166 * t202 - t167 * t203) * MDP(12) + (-t132 * t203 + t133 * t202) * MDP(13)) * pkin(3) + (-t206 * t260 + t259) * t209 + ((-t141 + t145) * MDP(12) - t148 * MDP(14) + (t138 - t234) * MDP(15) - t139 * MDP(16)) * t175 + (t253 * t235 + (MDP(10) * t253 - pkin(3) * t243) * t204) * qJD(2); (t141 * t178 + t142 * t175 + t182) * MDP(13) + t190 * MDP(16) + (t140 * t175 + (-qJD(5) - t138) * t178 + t215) * MDP(17) + ((t178 - t256) * MDP(14) + (t175 - t225) * MDP(16)) * qJD(3) + t261 * (-t175 ^ 2 - t174); t178 * t236 + (-t190 + (t175 + t225) * qJD(3)) * MDP(15) + (-t174 - t208) * MDP(16) + (-qJD(3) * t140 + t139 * t178 + t132) * MDP(17);];
tauc = t1;
