% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:40:01
% EndTime: 2019-12-05 17:40:05
% DurationCPUTime: 1.33s
% Computational Cost: add. (850->162), mult. (1995->225), div. (0->0), fcn. (1454->6), ass. (0->85)
t189 = sin(pkin(8));
t190 = cos(pkin(8));
t193 = sin(qJ(4));
t195 = cos(qJ(4));
t170 = -t193 * t189 + t195 * t190;
t169 = t195 * t189 + t193 * t190;
t167 = t169 * qJD(4);
t159 = qJD(1) * t167;
t223 = qJD(1) * t193;
t216 = t189 * t223;
t175 = qJD(4) * t216;
t222 = qJD(1) * t195;
t215 = t190 * t222;
t160 = qJD(4) * t215 - t175;
t166 = t215 - t216;
t192 = sin(qJ(5));
t194 = cos(qJ(5));
t221 = qJD(5) * t192;
t224 = qJD(1) * t169;
t227 = t194 * t224;
t125 = -qJD(5) * t227 - t194 * t159 - t192 * t160 - t166 * t221;
t205 = t194 * t166 - t192 * t224;
t126 = t205 * qJD(5) - t192 * t159 + t194 * t160;
t140 = t192 * t166 + t227;
t186 = qJD(4) + qJD(5);
t233 = t140 * t186;
t234 = t205 * t186;
t248 = t140 * t205 * MDP(18) + (-t126 + t234) * MDP(21) + (-t140 ^ 2 + t205 ^ 2) * MDP(19) + (t125 + t233) * MDP(20);
t247 = t170 * qJD(3);
t225 = t189 ^ 2 + t190 ^ 2;
t246 = t225 * qJD(3);
t202 = t247 * qJD(1);
t191 = -pkin(1) - qJ(3);
t240 = t191 * qJD(1);
t176 = qJD(2) + t240;
t212 = -pkin(6) * qJD(1) + t176;
t161 = t212 * t189;
t162 = t212 * t190;
t206 = -t195 * t161 - t193 * t162;
t128 = t159 * pkin(7) + t206 * qJD(4) - t202;
t136 = -pkin(7) * t224 - t206;
t188 = qJD(1) * qJ(2);
t182 = qJD(3) + t188;
t183 = t189 * pkin(3);
t173 = qJD(1) * t183 + t182;
t146 = pkin(4) * t224 + t173;
t245 = t146 * t140 + t136 * t221 + (-t136 * t186 - t128) * t192;
t244 = qJ(2) * MDP(6) + t189 * MDP(7) + t190 * MDP(8) + MDP(5);
t241 = -t193 * t161 + t195 * t162;
t239 = qJD(5) - t186;
t200 = t169 * qJD(3);
t127 = -t160 * pkin(7) - qJD(1) * t200 + t241 * qJD(4);
t238 = -t192 * t127 + t194 * t128 - t146 * t205;
t237 = t166 * pkin(4);
t236 = -pkin(6) + t191;
t228 = t194 * t136;
t180 = qJ(2) + t183;
t219 = t193 * qJD(4);
t218 = t195 * qJD(4);
t187 = qJD(1) * qJD(2);
t135 = -t166 * pkin(7) + t241;
t134 = qJD(4) * pkin(4) + t135;
t213 = -pkin(4) * t186 - t134;
t168 = -t189 * t219 + t190 * t218;
t209 = -t194 * t167 - t192 * t168;
t208 = qJD(1) * t225;
t204 = -t192 * t167 + t194 * t168;
t144 = t194 * t169 + t192 * t170;
t145 = -t192 * t169 + t194 * t170;
t171 = t236 * t189;
t172 = t236 * t190;
t203 = -t195 * t171 - t193 * t172;
t198 = -t171 * t219 + t172 * t218 - t200;
t197 = t203 * qJD(4) - t247;
t196 = qJD(1) ^ 2;
t154 = t168 * pkin(4) + qJD(2);
t152 = t169 * pkin(4) + t180;
t147 = t160 * pkin(4) + t187;
t138 = -t169 * pkin(7) - t203;
t137 = -t170 * pkin(7) - t193 * t171 + t195 * t172;
t133 = t167 * pkin(7) + t197;
t132 = -t168 * pkin(7) + t198;
t130 = t145 * qJD(5) + t204;
t129 = -t144 * qJD(5) + t209;
t1 = [0.2e1 * qJD(3) * MDP(9) * t208 + ((t182 + t188) * qJD(2) + (-t176 - t240) * t246) * MDP(10) + (-t159 * t170 - t166 * t167) * MDP(11) + (t159 * t169 - t170 * t160 - t166 * t168 + t167 * t224) * MDP(12) + (0.2e1 * t224 * qJD(2) + t180 * t160 + t173 * t168) * MDP(16) + (-t180 * t159 - t173 * t167 + (qJD(1) * t170 + t166) * qJD(2)) * MDP(17) + (t125 * t145 + t129 * t205) * MDP(18) + (-t125 * t144 - t145 * t126 - t129 * t140 - t130 * t205) * MDP(19) + (t152 * t126 + t146 * t130 + t154 * t140 + t147 * t144) * MDP(23) + (t152 * t125 + t146 * t129 + t147 * t145 + t154 * t205) * MDP(24) + (t129 * MDP(20) - t130 * MDP(21) + (-t192 * t132 + t194 * t133 + (-t137 * t192 - t138 * t194) * qJD(5)) * MDP(23) + (-t194 * t132 - t192 * t133 - (t137 * t194 - t138 * t192) * qJD(5)) * MDP(24)) * t186 + (-t167 * MDP(13) - t168 * MDP(14) + t197 * MDP(16) - t198 * MDP(17)) * qJD(4) + 0.2e1 * t244 * t187; (t209 * MDP(23) - t204 * MDP(24) + (-t144 * MDP(23) - t145 * MDP(24)) * qJD(5)) * t186 + (-t167 * MDP(16) - t168 * MDP(17)) * qJD(4) + ((-t182 - t246) * MDP(10) - t224 * MDP(16) - t166 * MDP(17) - t140 * MDP(23) - t205 * MDP(24)) * qJD(1) - t244 * t196; (t176 * t208 + t187) * MDP(10) - t175 * MDP(16) + (t126 + t234) * MDP(23) + (t125 - t233) * MDP(24) - t225 * MDP(9) * t196 + ((t166 + t215) * MDP(16) + (-t189 * t222 - t190 * t223 - t224) * MDP(17)) * qJD(4); t166 ^ 2 * MDP(12) + (t175 + (t166 - t215) * qJD(4)) * MDP(14) + (-t173 * t166 - t202) * MDP(16) + (-t140 * t237 - (-t192 * t135 - t228) * t186 + (t213 * t192 - t228) * qJD(5) + t238) * MDP(23) + (-t205 * t237 + (t213 * qJD(5) + t135 * t186 - t127) * t194 + t245) * MDP(24) + (t166 * MDP(11) + (qJD(3) + t173) * MDP(17) - MDP(12) * t224) * t224 + t248; (t239 * (-t192 * t134 - t228) + t238) * MDP(23) + ((-t239 * t134 - t127) * t194 + t245) * MDP(24) + t248;];
tauc = t1;
