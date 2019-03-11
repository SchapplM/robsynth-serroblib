% Calculate joint inertia matrix for
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR12_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:41:47
% EndTime: 2019-03-09 14:41:50
% DurationCPUTime: 1.23s
% Computational Cost: add. (1352->231), mult. (2898->319), div. (0->0), fcn. (3161->10), ass. (0->106)
t200 = sin(pkin(6));
t205 = sin(qJ(2));
t250 = t200 * t205;
t183 = pkin(8) * t250;
t201 = cos(pkin(6));
t208 = cos(qJ(2));
t259 = pkin(1) * t208;
t233 = -pkin(2) - t259;
t150 = pkin(3) * t250 + t183 + (-pkin(9) + t233) * t201;
t209 = -pkin(2) - pkin(9);
t228 = -qJ(3) * t205 - pkin(1);
t156 = (t209 * t208 + t228) * t200;
t204 = sin(qJ(4));
t207 = cos(qJ(4));
t134 = t207 * t150 - t156 * t204;
t135 = t150 * t204 + t156 * t207;
t249 = t200 * t208;
t165 = t201 * t204 + t207 * t249;
t166 = t201 * t207 - t204 * t249;
t272 = t166 * MDP(17) - t165 * MDP(18) + t134 * MDP(20) - t135 * MDP(21);
t202 = sin(qJ(6));
t206 = cos(qJ(6));
t247 = t202 * MDP(31) + t206 * MDP(32);
t270 = t202 * MDP(34) + t206 * MDP(35);
t269 = -MDP(34) * t206 + MDP(35) * t202;
t203 = sin(qJ(5));
t260 = cos(qJ(5));
t144 = t260 * t165 + t166 * t203;
t145 = -t203 * t165 + t260 * t166;
t268 = t145 * MDP(24) - t144 * MDP(25);
t176 = t203 * t204 - t260 * t207;
t177 = t203 * t207 + t260 * t204;
t267 = -t176 * MDP(27) - t177 * MDP(28);
t257 = -pkin(10) + t209;
t179 = t257 * t204;
t229 = t257 * t207;
t153 = t179 * t203 - t260 * t229;
t154 = t260 * t179 + t203 * t229;
t265 = -t176 * MDP(24) - t177 * MDP(25) - t153 * MDP(27) - t154 * MDP(28);
t264 = t176 ^ 2;
t174 = t177 ^ 2;
t263 = -2 * MDP(16);
t262 = 0.2e1 * MDP(34);
t261 = 0.2e1 * MDP(35);
t256 = pkin(2) * MDP(14);
t235 = pkin(4) * t250;
t130 = -pkin(10) * t166 + t134 + t235;
t131 = -pkin(10) * t165 + t135;
t125 = t260 * t130 - t203 * t131;
t123 = -pkin(5) * t250 - t125;
t255 = t123 * t206;
t254 = t144 * t202;
t253 = t144 * t206;
t252 = t153 * t206;
t251 = t176 * t202;
t248 = t201 * t208;
t187 = t204 * pkin(4) + qJ(3);
t168 = t201 * t205 * pkin(1) + pkin(8) * t249;
t246 = MDP(23) * t145;
t138 = t145 * t206 + t202 * t250;
t245 = MDP(31) * t138;
t137 = t145 * t202 - t206 * t250;
t244 = MDP(32) * t137;
t243 = MDP(33) * t144;
t240 = t138 * MDP(29);
t238 = t206 * MDP(29);
t236 = t207 * MDP(15);
t234 = t260 * pkin(4);
t191 = t201 * qJ(3);
t157 = -t191 - t168;
t232 = t260 * t131;
t231 = t202 * t206 * MDP(30);
t197 = t202 ^ 2;
t230 = t197 * MDP(29) + MDP(26) + 0.2e1 * t231;
t155 = pkin(3) * t249 - t157;
t227 = pkin(5) * t176 - pkin(11) * t177;
t189 = pkin(4) * t203 + pkin(11);
t190 = -t234 - pkin(5);
t226 = -t176 * t190 - t177 * t189;
t225 = -t168 * MDP(10) + (pkin(1) * t248 - t183) * MDP(9);
t222 = MDP(20) * t207 - MDP(21) * t204;
t126 = t203 * t130 + t232;
t221 = t125 * MDP(27) - t126 * MDP(28);
t220 = t144 * MDP(27) + t145 * MDP(28);
t219 = -t206 * MDP(31) + t202 * MDP(32);
t139 = pkin(4) * t165 + t155;
t216 = MDP(23) + t219;
t215 = t269 * t176 + t267;
t214 = (-t137 * t202 + t138 * t206) * MDP(30) + t202 * t240 + MDP(26) * t250 + t268 + t247 * t144;
t213 = (t260 * MDP(27) - t203 * MDP(28)) * pkin(4);
t199 = t206 ^ 2;
t212 = (t197 - t199) * t176 * MDP(30) - t238 * t251 + t265 + t247 * t177;
t211 = (MDP(20) * t209 + MDP(17)) * t207 + (-MDP(21) * t209 - MDP(18)) * t204;
t124 = pkin(11) * t250 + t126;
t127 = pkin(5) * t144 - pkin(11) * t145 + t139;
t120 = -t124 * t202 + t127 * t206;
t121 = t124 * t206 + t127 * t202;
t210 = MDP(34) * t120 - MDP(35) * t121 + t243 - t244 + t245;
t160 = t233 * t201 + t183;
t158 = (-pkin(2) * t208 + t228) * t200;
t149 = t153 * t202;
t148 = pkin(5) * t177 + pkin(11) * t176 + t187;
t133 = t148 * t202 + t154 * t206;
t132 = t148 * t206 - t154 * t202;
t122 = t123 * t202;
t1 = [t145 ^ 2 * MDP(22) + (t157 ^ 2 + t158 ^ 2 + t160 ^ 2) * MDP(14) + t201 ^ 2 * MDP(8) + MDP(1) + (t166 * MDP(15) + t165 * t263) * t166 + (-0.2e1 * t137 * MDP(30) + t240) * t138 + (t243 - 0.2e1 * t244 + 0.2e1 * t245 - 0.2e1 * t246) * t144 + (t120 * t144 + t123 * t137) * t262 + (-t121 * t144 + t123 * t138) * t261 + 0.2e1 * (t160 * MDP(12) - t157 * MDP(13) + t225) * t201 + 0.2e1 * (t165 * MDP(20) + t166 * MDP(21)) * t155 + 0.2e1 * t220 * t139 + ((0.2e1 * t208 * MDP(5) + (MDP(19) + MDP(26) + MDP(4)) * t205) * t205 + 0.2e1 * (-MDP(10) * t205 + MDP(9) * t208) * pkin(1)) * t200 ^ 2 + 0.2e1 * (MDP(7) * t248 + (-t157 * MDP(11) + t158 * MDP(12)) * t208 + (t160 * MDP(11) - t158 * MDP(13) + MDP(6) * t201 + t221 + t268 + t272) * t205) * t200; (-pkin(2) * t160 - qJ(3) * t157) * MDP(14) + t166 * t236 + (qJ(3) * t165 + t155 * t204) * MDP(20) + (qJ(3) * t166 + t155 * t207) * MDP(21) + (0.2e1 * t191 + t168) * MDP(13) + (t132 * t144 + t137 * t153) * MDP(34) + (-t133 * t144 + t138 * t153) * MDP(35) + (-t165 * t207 - t166 * t204) * MDP(16) + t183 * MDP(12) + (MDP(8) + (-0.2e1 * pkin(2) - t259) * MDP(12)) * t201 + t220 * t187 + (MDP(27) * t139 + t210 - t246) * t177 + (-t138 * t238 - t139 * MDP(28) + (t137 * t206 + t138 * t202) * MDP(30) - t145 * MDP(22) - t270 * t123 + t216 * t144) * t176 + ((qJ(3) * MDP(11) + MDP(7)) * t208 + (-pkin(2) * MDP(11) + MDP(6) + t211 + t265) * t205) * t200 + t225; 0.2e1 * t187 * t177 * MDP(27) + t174 * MDP(33) + MDP(8) + (t204 * t263 + t236) * t207 + (-0.2e1 * MDP(12) + t256) * pkin(2) + 0.2e1 * (-t187 * MDP(28) + t177 * t216) * t176 + (t132 * t177 - t153 * t251) * t262 + (-t133 * t177 - t176 * t252) * t261 + (MDP(29) * t199 + MDP(22) - 0.2e1 * t231) * t264 + (MDP(14) * qJ(3) + 0.2e1 * t204 * MDP(20) + 0.2e1 * t207 * MDP(21) + 0.2e1 * MDP(13)) * qJ(3); t201 * MDP(12) + t160 * MDP(14) + (t137 * t176 - t177 * t254) * MDP(34) + (t138 * t176 - t177 * t253) * MDP(35) + (MDP(11) + t222 + t267) * t250; MDP(12) - t256 + t270 * (-t174 - t264); MDP(14); MDP(19) * t250 + (t234 * t250 + t125) * MDP(27) + (-t232 + (-t130 - t235) * t203) * MDP(28) + (t137 * t190 - t189 * t254 - t255) * MDP(34) + (t138 * t190 - t189 * t253 + t122) * MDP(35) + t214 + t272; (t226 * t202 - t252) * MDP(34) + (t226 * t206 + t149) * MDP(35) + t211 + t212; t215 + t222; 0.2e1 * t190 * t269 + MDP(19) + 0.2e1 * t213 + t230; (-pkin(5) * t137 - pkin(11) * t254 - t255) * MDP(34) + (-pkin(5) * t138 - pkin(11) * t253 + t122) * MDP(35) + t214 + t221; (t227 * t202 - t252) * MDP(34) + (t227 * t206 + t149) * MDP(35) + t212; t215; t213 + t230 - t269 * (pkin(5) - t190); -0.2e1 * pkin(5) * t269 + t230; t210; MDP(33) * t177 + t132 * MDP(34) - t133 * MDP(35) + t219 * t176; -t270 * t177; -t189 * t270 + t247; -pkin(11) * t270 + t247; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
