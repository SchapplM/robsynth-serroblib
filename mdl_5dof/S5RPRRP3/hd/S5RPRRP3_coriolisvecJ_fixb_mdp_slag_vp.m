% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:03:57
% EndTime: 2019-12-05 18:04:01
% DurationCPUTime: 0.99s
% Computational Cost: add. (1033->150), mult. (2447->212), div. (0->0), fcn. (1560->6), ass. (0->89)
t192 = sin(qJ(3));
t193 = cos(qJ(3));
t181 = sin(pkin(8)) * pkin(1) + pkin(6);
t241 = pkin(7) + t181;
t212 = t241 * qJD(1);
t154 = qJD(2) * t192 + t193 * t212;
t191 = sin(qJ(4));
t149 = t191 * t154;
t153 = t193 * qJD(2) - t212 * t192;
t240 = qJD(3) * pkin(3);
t152 = t153 + t240;
t242 = cos(qJ(4));
t210 = t242 * t152 - t149;
t171 = t191 * t193 + t192 * t242;
t227 = qJD(1) * t171;
t237 = t227 * qJ(5);
t247 = t237 - t210;
t186 = qJD(3) + qJD(4);
t246 = MDP(5) * t193;
t245 = MDP(6) * (t192 ^ 2 - t193 ^ 2);
t244 = -t192 * MDP(10) - t193 * MDP(11);
t243 = t227 ^ 2;
t218 = t242 * t193;
t208 = qJD(1) * t218;
t225 = qJD(1) * t192;
t217 = t191 * t225;
t163 = -t208 + t217;
t182 = -cos(pkin(8)) * pkin(1) - pkin(2);
t172 = -pkin(3) * t193 + t182;
t167 = t172 * qJD(1);
t141 = pkin(4) * t163 + qJD(5) + t167;
t239 = t141 * t227;
t238 = t163 * qJ(5);
t236 = t167 * t227;
t235 = t191 * t192;
t194 = qJD(3) ^ 2;
t234 = t192 * t194;
t233 = t193 * t194;
t127 = pkin(4) * t186 - t247;
t232 = t127 + t247;
t145 = t186 * t171;
t140 = t145 * qJD(1);
t207 = t186 * t235;
t216 = qJD(4) * t242;
t144 = -qJD(3) * t218 - t193 * t216 + t207;
t231 = -t171 * t140 + t144 * t163;
t230 = t242 * t153 - t149;
t229 = t186 * t208;
t174 = qJD(1) * t182;
t223 = qJD(4) * t191;
t220 = qJD(1) * qJD(3);
t219 = t192 * t240;
t151 = t242 * t154;
t215 = t192 * t220;
t180 = pkin(3) * t215;
t214 = pkin(4) * t140 + t180;
t213 = qJD(3) * t241;
t147 = t153 * qJD(3);
t148 = t154 * qJD(3);
t211 = -t191 * t147 - t242 * t148;
t209 = -t153 * t191 - t151;
t139 = qJD(1) * t207 - t229;
t170 = -t218 + t235;
t205 = -t139 * t170 + t145 * t227;
t204 = 0.2e1 * qJD(3) * t174;
t203 = -t191 * t152 - t151;
t168 = t241 * t192;
t169 = t241 * t193;
t202 = t191 * t168 - t169 * t242;
t162 = t163 ^ 2;
t201 = t227 * t163 * MDP(12) + (-t162 + t243) * MDP(13) + (t229 + (t163 - t217) * t186) * MDP(14);
t160 = t192 * t213;
t161 = t193 * t213;
t200 = -t242 * t160 - t191 * t161 - t168 * t216 - t169 * t223;
t199 = qJD(4) * t203 + t211;
t198 = qJD(4) * t202 + t191 * t160 - t242 * t161;
t197 = t147 * t242 - t191 * t148 + t152 * t216 - t154 * t223;
t196 = t167 * t163 - t197;
t184 = pkin(3) * t242 + pkin(4);
t135 = -t170 * qJ(5) - t202;
t134 = -t171 * qJ(5) - t168 * t242 - t191 * t169;
t131 = -t237 + t230;
t130 = t209 + t238;
t129 = -t203 - t238;
t126 = t144 * qJ(5) - t171 * qJD(5) + t198;
t125 = -qJ(5) * t145 - qJD(5) * t170 + t200;
t124 = t139 * qJ(5) - qJD(5) * t227 + t199;
t123 = -t140 * qJ(5) - t163 * qJD(5) + t197;
t1 = [0.2e1 * t215 * t246 - 0.2e1 * t220 * t245 + MDP(7) * t233 - MDP(8) * t234 + (-t181 * t233 + t192 * t204) * MDP(10) + (t181 * t234 + t193 * t204) * MDP(11) + (-t139 * t171 - t144 * t227) * MDP(12) + (-t205 + t231) * MDP(13) + (t172 * t140 + t167 * t145 + t163 * t219 + t180 * t170) * MDP(17) + (-t172 * t139 - t167 * t144 + 0.2e1 * t227 * t219) * MDP(18) + (-t123 * t170 - t124 * t171 - t125 * t163 - t126 * t227 + t127 * t144 - t129 * t145 + t134 * t139 - t135 * t140) * MDP(19) + (t123 * t135 + t129 * t125 + t124 * t134 + t127 * t126 + t214 * (pkin(4) * t170 + t172) + t141 * (pkin(4) * t145 + t219)) * MDP(20) + (-t144 * MDP(14) - t145 * MDP(15) + t198 * MDP(17) - t200 * MDP(18)) * t186; (t205 + t231) * MDP(19) + (t123 * t171 - t124 * t170 - t127 * t145 - t129 * t144) * MDP(20) + t244 * t194 + (-t145 * MDP(17) + t144 * MDP(18)) * t186; (-pkin(3) * t163 * t225 - t236 - t209 * t186 + (-t151 + (-pkin(3) * t186 - t152) * t191) * qJD(4) + t211) * MDP(17) + (t230 * t186 + (-t186 * t216 - t225 * t227) * pkin(3) + t196) * MDP(18) + (t184 * t139 + (t129 + t130) * t227 + (-t127 + t131) * t163 + (-t140 * t191 + (-t163 * t242 + t191 * t227) * qJD(4)) * pkin(3)) * MDP(19) + (-pkin(4) * t239 + t124 * t184 - t127 * t130 - t129 * t131 + (-t141 * t225 + t123 * t191 + (-t127 * t191 + t129 * t242) * qJD(4)) * pkin(3)) * MDP(20) + t201 + (t244 * t174 + (-t192 * t246 + t245) * qJD(1)) * qJD(1); (-t186 * t203 + t199 - t236) * MDP(17) + (t186 * t210 + t196) * MDP(18) + (pkin(4) * t139 - t163 * t232) * MDP(19) + (t232 * t129 + (t124 - t239) * pkin(4)) * MDP(20) + t201; (-t162 - t243) * MDP(19) + (t127 * t227 + t129 * t163 + t214) * MDP(20);];
tauc = t1;
