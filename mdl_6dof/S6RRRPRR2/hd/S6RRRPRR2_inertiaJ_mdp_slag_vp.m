% Calculate joint inertia matrix for
% S6RRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR2_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:08:42
% EndTime: 2019-03-09 18:08:43
% DurationCPUTime: 0.81s
% Computational Cost: add. (1444->172), mult. (2688->247), div. (0->0), fcn. (3185->10), ass. (0->94)
t198 = sin(qJ(5));
t202 = cos(qJ(5));
t197 = sin(qJ(6));
t201 = cos(qJ(6));
t174 = t197 * t198 - t201 * t202;
t175 = t197 * t202 + t201 * t198;
t229 = t175 * MDP(29) - t174 * MDP(30);
t222 = t198 * MDP(22) + t202 * MDP(23) + t229;
t215 = t202 * MDP(25) - t198 * MDP(26);
t245 = t174 * MDP(32) + t175 * MDP(33);
t249 = t215 - t245;
t199 = sin(qJ(3));
t200 = sin(qJ(2));
t203 = cos(qJ(3));
t204 = cos(qJ(2));
t176 = t199 * t204 + t203 * t200;
t195 = sin(pkin(11));
t196 = cos(pkin(11));
t216 = t199 * t200 - t203 * t204;
t152 = t196 * t176 - t195 * t216;
t130 = t175 * t152;
t131 = t174 * t152;
t248 = -t131 * MDP(29) - t130 * MDP(30);
t187 = t203 * pkin(2) + pkin(3);
t241 = pkin(2) * t199;
t163 = t195 * t187 + t196 * t241;
t161 = pkin(9) + t163;
t156 = (-pkin(10) - t161) * t198;
t192 = t202 * pkin(10);
t157 = t202 * t161 + t192;
t135 = t201 * t156 - t197 * t157;
t136 = t197 * t156 + t201 * t157;
t247 = t135 * MDP(32) - t136 * MDP(33);
t184 = t195 * pkin(3) + pkin(9);
t172 = (-pkin(10) - t184) * t198;
t173 = t202 * t184 + t192;
t146 = t201 * t172 - t197 * t173;
t147 = t197 * t172 + t201 * t173;
t246 = t146 * MDP(32) - t147 * MDP(33);
t188 = -t204 * pkin(2) - pkin(1);
t244 = 0.2e1 * t188;
t243 = -2 * MDP(28);
t242 = pkin(7) + pkin(8);
t151 = t195 * t176 + t196 * t216;
t240 = t151 * pkin(5);
t239 = t202 * pkin(5);
t179 = t242 * t200;
t180 = t242 * t204;
t217 = t199 * t179 - t203 * t180;
t141 = -t216 * qJ(4) - t217;
t218 = -t203 * t179 - t199 * t180;
t207 = -t176 * qJ(4) + t218;
t125 = t195 * t141 - t196 * t207;
t238 = t125 * t202;
t237 = t152 * t198;
t236 = t152 * t202;
t235 = t198 * t202;
t159 = pkin(3) * t216 + t188;
t124 = t151 * pkin(4) - t152 * pkin(9) + t159;
t127 = t196 * t141 + t195 * t207;
t233 = t202 * t127;
t114 = t233 + (-pkin(10) * t152 + t124) * t198;
t234 = t201 * t114;
t228 = MDP(27) * t131;
t225 = MDP(24) + MDP(31);
t224 = t151 * MDP(31) + t248;
t185 = -t196 * pkin(3) - pkin(4);
t223 = MDP(21) * t235;
t116 = t202 * t124 - t198 * t127;
t113 = -pkin(10) * t236 + t116 + t240;
t110 = t201 * t113 - t197 * t114;
t162 = t196 * t187 - t195 * t241;
t160 = -pkin(4) - t162;
t193 = t198 ^ 2;
t221 = t193 * MDP(20) + MDP(15) + 0.2e1 * t223 + (MDP(27) * t175 + t174 * t243) * t175;
t220 = -t151 * t161 + t152 * t160;
t219 = -t151 * t184 + t152 * t185;
t214 = -MDP(25) * t198 - MDP(26) * t202;
t212 = (t203 * MDP(16) - t199 * MDP(17)) * pkin(2);
t211 = (MDP(32) * t201 - MDP(33) * t197) * pkin(5);
t210 = (t202 * MDP(22) - t198 * MDP(23)) * t152;
t209 = 0.2e1 * t215;
t208 = 0.2e1 * t245;
t194 = t202 ^ 2;
t206 = (-t175 * t130 + t131 * t174) * MDP(28) - t175 * t228 + t218 * MDP(16) + t217 * MDP(17) - t216 * MDP(14) + t176 * MDP(13) + ((-t193 + t194) * MDP(21) + MDP(20) * t235) * t152 + t222 * t151;
t177 = t185 - t239;
t158 = t160 - t239;
t123 = t125 * t198;
t120 = pkin(5) * t237 + t125;
t119 = t120 * t175;
t118 = t120 * t174;
t117 = t198 * t124 + t233;
t111 = t197 * t113 + t234;
t1 = [(t125 ^ 2 + t127 ^ 2 + t159 ^ 2) * MDP(19) + MDP(1) + t216 * MDP(16) * t244 + (t194 * MDP(20) - 0.2e1 * t223) * t152 ^ 2 + t225 * t151 ^ 2 - (t130 * t243 - t228) * t131 + (MDP(11) * t176 - 0.2e1 * MDP(12) * t216 + MDP(17) * t244) * t176 + 0.2e1 * (t110 * t151 + t120 * t130) * MDP(32) + 0.2e1 * (-t111 * t151 - t120 * t131) * MDP(33) + 0.2e1 * (t125 * t152 - t127 * t151) * MDP(18) + 0.2e1 * (t116 * t151 + t125 * t237) * MDP(25) + 0.2e1 * (-t117 * t151 + t125 * t236) * MDP(26) + (MDP(4) * t200 + 0.2e1 * t204 * MDP(5)) * t200 + 0.2e1 * (-t200 * MDP(10) + t204 * MDP(9)) * pkin(1) + 0.2e1 * (t248 + t210) * t151; (t220 * t202 + t123) * MDP(26) + (t158 * t130 + t135 * t151 + t118) * MDP(32) + (-t158 * t131 - t136 * t151 + t119) * MDP(33) + t206 + (t220 * t198 - t238) * MDP(25) + (-t204 * MDP(10) - t200 * MDP(9)) * pkin(7) + (-t163 * t151 - t162 * t152) * MDP(18) + (-t125 * t162 + t127 * t163) * MDP(19) + t200 * MDP(6) + t204 * MDP(7); MDP(8) + (t162 ^ 2 + t163 ^ 2) * MDP(19) - t160 * t209 + t158 * t208 + 0.2e1 * t212 + t221; t206 + (t219 * t202 + t123) * MDP(26) + (t219 * t198 - t238) * MDP(25) + ((-t151 * t195 - t152 * t196) * MDP(18) + (-t125 * t196 + t127 * t195) * MDP(19)) * pkin(3) + (t177 * t130 + t146 * t151 + t118) * MDP(32) + (-t177 * t131 - t147 * t151 + t119) * MDP(33); (t162 * t196 + t163 * t195) * MDP(19) * pkin(3) + t212 + t221 + t245 * (t158 + t177) - t215 * (t160 + t185); (t195 ^ 2 + t196 ^ 2) * MDP(19) * pkin(3) ^ 2 - t185 * t209 + t177 * t208 + t221; t159 * MDP(19) + t249 * t151; 0; 0; MDP(19); t151 * MDP(24) + t116 * MDP(25) - t117 * MDP(26) + (t201 * t240 + t110) * MDP(32) + (-t234 + (-t113 - t240) * t197) * MDP(33) + t210 + t224; t161 * t214 + t222 + t247; t184 * t214 + t222 + t246; t249; 0.2e1 * t211 + t225; t110 * MDP(32) - t111 * MDP(33) + t224; t229 + t247; t229 + t246; -t245; MDP(31) + t211; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
