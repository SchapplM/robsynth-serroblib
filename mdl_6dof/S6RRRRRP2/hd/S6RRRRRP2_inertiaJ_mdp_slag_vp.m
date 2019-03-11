% Calculate joint inertia matrix for
% S6RRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRRP2_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:03:14
% EndTime: 2019-03-10 01:03:18
% DurationCPUTime: 1.07s
% Computational Cost: add. (1882->212), mult. (3348->272), div. (0->0), fcn. (3753->8), ass. (0->104)
t195 = cos(qJ(3));
t178 = pkin(2) * t195 + pkin(3);
t190 = sin(qJ(4));
t194 = cos(qJ(4));
t191 = sin(qJ(3));
t257 = pkin(2) * t191;
t158 = -t178 * t190 - t194 * t257;
t156 = pkin(10) - t158;
t189 = sin(qJ(5));
t187 = t189 ^ 2;
t193 = cos(qJ(5));
t188 = t193 ^ 2;
t241 = t187 + t188;
t246 = t241 * t156;
t176 = pkin(3) * t190 + pkin(10);
t266 = t241 * t176;
t265 = t189 * MDP(27) + t193 * MDP(28);
t262 = t193 * MDP(30);
t263 = t189 * MDP(31);
t201 = 0.2e1 * t262 - 0.2e1 * t263;
t196 = cos(qJ(2));
t179 = -t196 * pkin(2) - pkin(1);
t192 = sin(qJ(2));
t209 = t191 * t192 - t195 * t196;
t153 = t209 * pkin(3) + t179;
t261 = 0.2e1 * t153;
t260 = 0.2e1 * t179;
t259 = 2 * MDP(33);
t258 = pkin(7) + pkin(8);
t256 = pkin(3) * t194;
t163 = t191 * t196 + t192 * t195;
t145 = t163 * t190 + t194 * t209;
t255 = pkin(5) * t145;
t254 = pkin(10) * t145;
t253 = pkin(4) * MDP(31);
t252 = qJ(6) * t145;
t251 = t145 * t156;
t250 = t145 * t176;
t249 = t189 * t193;
t146 = t194 * t163 - t190 * t209;
t133 = t145 * pkin(4) - t146 * pkin(10) + t153;
t167 = t258 * t192;
t168 = t258 * t196;
t221 = -t195 * t167 - t168 * t191;
t138 = -pkin(9) * t163 + t221;
t210 = t191 * t167 - t195 * t168;
t139 = -t209 * pkin(9) - t210;
t135 = t138 * t190 + t139 * t194;
t125 = t189 * t133 + t193 * t135;
t217 = -pkin(5) * t193 - qJ(6) * t189;
t166 = -pkin(4) + t217;
t243 = -t194 * t178 + t190 * t257;
t149 = t166 + t243;
t162 = t166 - t256;
t248 = -t149 - t162;
t247 = -t149 - t166;
t245 = -t162 - t166;
t242 = t241 * pkin(10);
t240 = MDP(24) * t190;
t239 = MDP(34) * t189;
t238 = MDP(35) * t162;
t237 = MDP(35) * t189;
t122 = t125 + t252;
t236 = t122 * MDP(35);
t222 = -t193 * t133 + t135 * t189;
t123 = t222 - t255;
t235 = t123 * MDP(35);
t134 = -t138 * t194 + t139 * t190;
t216 = pkin(5) * t189 - qJ(6) * t193;
t126 = t216 * t146 + t134;
t234 = t126 * MDP(34);
t233 = t145 * MDP(29);
t232 = t149 * MDP(35);
t231 = t243 * MDP(23);
t230 = t158 * MDP(24);
t229 = t166 * MDP(35);
t228 = t195 * MDP(16);
t227 = -t216 * MDP(33) + t265;
t226 = MDP(26) * t249;
t225 = t187 * MDP(25) + MDP(22) + 0.2e1 * t226;
t223 = -MDP(35) * pkin(5) - MDP(32);
t220 = t241 * MDP(35);
t219 = MDP(15) + t225;
t218 = -pkin(4) * t146 - t254;
t215 = -t146 * t166 + t254;
t214 = -t146 * t149 + t251;
t155 = -pkin(4) + t243;
t213 = t146 * t155 - t251;
t212 = -t146 * t162 + t250;
t177 = -pkin(4) - t256;
t211 = t146 * t177 - t250;
t208 = t193 * MDP(27) - t189 * MDP(28);
t207 = -t222 * MDP(30) - t125 * MDP(31);
t205 = -t134 * MDP(30) - t126 * MDP(32);
t204 = -0.2e1 * MDP(32) * t193 - 0.2e1 * t239;
t203 = pkin(4) * t262 + t225;
t202 = (MDP(23) * t194 - t240) * pkin(3);
t200 = (t122 * t193 + t123 * t189) * MDP(33) - t135 * MDP(24) + (t263 - MDP(23)) * t134 + ((-t187 + t188) * MDP(26) + MDP(25) * t249 + MDP(20)) * t146 + (-MDP(21) + t265) * t145;
t199 = (MDP(35) * qJ(6) - MDP(31) + MDP(34)) * t193 + (-MDP(30) + t223) * t189;
t198 = t163 * MDP(13) - t209 * MDP(14) + t221 * MDP(16) + t210 * MDP(17) + t200;
t181 = t189 * MDP(33);
t172 = t177 * t189;
t152 = t155 * t189;
t1 = [(t122 ^ 2 + t123 ^ 2 + t126 ^ 2) * MDP(35) + t209 * MDP(16) * t260 + MDP(1) + (MDP(11) * t163 - 0.2e1 * t209 * MDP(12) + MDP(17) * t260) * t163 + (MDP(23) * t261 + t233) * t145 + 0.2e1 * (-t123 * MDP(32) + t122 * MDP(34) + t207) * t145 + (MDP(4) * t192 + 0.2e1 * t196 * MDP(5)) * t192 + 0.2e1 * (-t192 * MDP(10) + t196 * MDP(9)) * pkin(1) + 0.2e1 * ((-t122 * t189 + t123 * t193) * MDP(33) + (t189 * MDP(30) + t193 * MDP(31)) * t134 + (t189 * MDP(32) - t193 * MDP(34)) * t126 + (-MDP(19) + t208) * t145) * t146 + (MDP(24) * t261 + (t188 * MDP(25) + MDP(18) - 0.2e1 * t226) * t146) * t146; (-t196 * MDP(10) - t192 * MDP(9)) * pkin(7) + t192 * MDP(6) + t196 * MDP(7) + t198 + t126 * t232 + (t213 * MDP(30) - t214 * MDP(32) + t156 * t235 - t234) * t189 + (t213 * MDP(31) + t214 * MDP(34) + t156 * t236 + t205) * t193; MDP(8) - t155 * t201 + t156 ^ 2 * t220 + (t204 + t232) * t149 + 0.2e1 * (-t191 * MDP(17) + t228) * pkin(2) - 0.2e1 * t231 + 0.2e1 * t230 + t246 * t259 + t219; (t211 * MDP(31) + t212 * MDP(34) + t176 * t236 + t205) * t193 + (t211 * MDP(30) - t212 * MDP(32) + t176 * t235 - t234) * t189 + t198 + t126 * t238; (-t243 + t256) * MDP(23) + (t172 + t152) * MDP(31) + (t266 + t246) * MDP(33) + (t149 * t162 + t156 * t266) * MDP(35) + (-pkin(3) - t178) * t240 + t248 * t239 + (t228 + (-MDP(24) * t194 - MDP(17)) * t191) * pkin(2) + ((-t155 - t177) * MDP(30) + t248 * MDP(32)) * t193 + t219; t266 * t259 - t177 * t201 + t176 ^ 2 * t220 + (t204 + t238) * t162 + 0.2e1 * t202 + t219; t126 * t229 + (t218 * MDP(31) + t215 * MDP(34) + pkin(10) * t236 + t205) * t193 + (t218 * MDP(30) - t215 * MDP(32) + pkin(10) * t235 - t234) * t189 + t200; -t231 + t230 + t152 * MDP(31) + (t242 + t246) * MDP(33) + (pkin(10) * t246 + t149 * t166) * MDP(35) + (-t155 * MDP(30) + t247 * MDP(32)) * t193 + (t247 * MDP(34) - t253) * t189 + t203; t172 * MDP(31) + (t242 + t266) * MDP(33) + (pkin(10) * t266 + t162 * t166) * MDP(35) + (-t177 * MDP(30) + t245 * MDP(32)) * t193 + (t245 * MDP(34) - t253) * t189 + t202 + t203; t242 * t259 + pkin(10) ^ 2 * t220 + (t204 + t229) * t166 + pkin(4) * t201 + t225; t233 + (-t222 + 0.2e1 * t255) * MDP(32) + (t125 + 0.2e1 * t252) * MDP(34) + (-pkin(5) * t123 + qJ(6) * t122) * MDP(35) + (t217 * MDP(33) + t208) * t146 + t207; t199 * t156 + t227; t199 * t176 + t227; t199 * pkin(10) + t227; MDP(29) + 0.2e1 * pkin(5) * MDP(32) + 0.2e1 * qJ(6) * MDP(34) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(35); t146 * t193 * MDP(33) - t145 * MDP(32) + t235; t156 * t237 + t181; t176 * t237 + t181; pkin(10) * t237 + t181; t223; MDP(35);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
