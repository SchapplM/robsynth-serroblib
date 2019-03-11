% Calculate joint inertia matrix for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR6_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:23:14
% EndTime: 2019-03-09 22:23:17
% DurationCPUTime: 1.14s
% Computational Cost: add. (2004->209), mult. (4001->305), div. (0->0), fcn. (4443->10), ass. (0->103)
t201 = sin(qJ(3));
t205 = cos(qJ(3));
t242 = pkin(8) + pkin(9);
t185 = t242 * t201;
t186 = t242 * t205;
t200 = sin(qJ(4));
t204 = cos(qJ(4));
t160 = -t204 * t185 - t200 * t186;
t161 = -t200 * t185 + t204 * t186;
t178 = t200 * t201 - t204 * t205;
t179 = t200 * t205 + t204 * t201;
t149 = -t179 * qJ(5) + t160;
t150 = -t178 * qJ(5) + t161;
t197 = sin(pkin(11));
t198 = cos(pkin(11));
t130 = t198 * t149 - t197 * t150;
t154 = -t197 * t178 + t198 * t179;
t120 = -t154 * pkin(10) + t130;
t131 = t197 * t149 + t198 * t150;
t153 = -t198 * t178 - t197 * t179;
t121 = t153 * pkin(10) + t131;
t199 = sin(qJ(6));
t203 = cos(qJ(6));
t136 = -t203 * t153 + t199 * t154;
t137 = t199 * t153 + t203 * t154;
t217 = t137 * MDP(29) - t136 * MDP(30) + (t203 * t120 - t199 * t121) * MDP(32) - (t199 * t120 + t203 * t121) * MDP(33);
t209 = t179 * MDP(20) - t178 * MDP(21) + t160 * MDP(23) - t161 * MDP(24) + t217;
t254 = t209 - (t201 * MDP(16) + t205 * MDP(17)) * pkin(8) + t201 * MDP(13) + t205 * MDP(14);
t210 = (t204 * MDP(23) - t200 * MDP(24)) * pkin(3);
t202 = sin(qJ(2));
t169 = t179 * t202;
t170 = t178 * t202;
t142 = -t198 * t169 + t197 * t170;
t143 = -t197 * t169 - t198 * t170;
t124 = -t203 * t142 + t199 * t143;
t125 = t199 * t142 + t203 * t143;
t233 = t125 * MDP(29) - t124 * MDP(30);
t253 = -t170 * MDP(20) - t169 * MDP(21);
t206 = cos(qJ(2));
t182 = -t206 * pkin(2) - t202 * pkin(8) - pkin(1);
t177 = t205 * t182;
t236 = pkin(9) * t202;
t239 = pkin(7) * t201;
t155 = -t205 * t236 + t177 + (-pkin(3) - t239) * t206;
t237 = pkin(7) * t206;
t220 = t205 * t237;
t157 = t220 + (t182 - t236) * t201;
t138 = t204 * t155 - t200 * t157;
t139 = t200 * t155 + t204 * t157;
t128 = -t206 * pkin(4) + t170 * qJ(5) + t138;
t132 = -t169 * qJ(5) + t139;
t118 = t198 * t128 - t197 * t132;
t116 = -t206 * pkin(5) - t143 * pkin(10) + t118;
t119 = t197 * t128 + t198 * t132;
t117 = t142 * pkin(10) + t119;
t109 = t203 * t116 - t199 * t117;
t110 = t199 * t116 + t203 * t117;
t251 = t109 * MDP(32) - t110 * MDP(33) + t233;
t252 = t138 * MDP(23) - t139 * MDP(24) + t251 + t253;
t248 = -2 * MDP(19);
t247 = 0.2e1 * MDP(23);
t246 = 0.2e1 * MDP(24);
t245 = 2 * MDP(25);
t244 = -2 * MDP(28);
t243 = 0.2e1 * MDP(33);
t241 = pkin(3) * t200;
t240 = pkin(4) * t197;
t238 = pkin(7) * t205;
t188 = t198 * pkin(4) + pkin(5);
t235 = t199 * t188;
t234 = t201 * t205;
t181 = (pkin(3) * t201 + pkin(7)) * t202;
t230 = t125 * MDP(27);
t229 = t136 * MDP(32);
t190 = t204 * pkin(3) + pkin(4);
t171 = t198 * t190 - t197 * t241;
t168 = pkin(5) + t171;
t173 = t197 * t190 + t198 * t241;
t147 = t203 * t168 - t199 * t173;
t226 = t147 * MDP(32);
t148 = t199 * t168 + t203 * t173;
t225 = t148 * MDP(33);
t224 = t170 * MDP(18);
t184 = t203 * t188;
t223 = (-t199 * t240 + t184) * MDP(32);
t222 = (t203 * t240 + t235) * MDP(33);
t221 = MDP(22) + MDP(31);
t191 = -t205 * pkin(3) - pkin(2);
t219 = MDP(12) * t234;
t218 = MDP(15) + t221;
t156 = t169 * pkin(4) + t181;
t165 = t178 * pkin(4) + t191;
t216 = t205 * MDP(13) - t201 * MDP(14);
t212 = MDP(31) - t225 + t226;
t211 = MDP(31) - t222 + t223;
t195 = t205 ^ 2;
t194 = t202 ^ 2;
t193 = t201 ^ 2;
t164 = t201 * t182 + t220;
t163 = -t201 * t237 + t177;
t140 = -t153 * pkin(5) + t165;
t133 = -t142 * pkin(5) + t156;
t1 = [-0.2e1 * pkin(1) * t202 * MDP(10) + MDP(1) + (t118 ^ 2 + t119 ^ 2 + t156 ^ 2) * MDP(26) - (t169 * t248 - t224) * t170 + (t124 * t244 + t230) * t125 + t218 * t206 ^ 2 + (t195 * MDP(11) + MDP(4) - 0.2e1 * t219) * t194 + 0.2e1 * (pkin(1) * MDP(9) + (MDP(5) - t216) * t202 - t253 - t233) * t206 + 0.2e1 * (t164 * t206 + t194 * t238) * MDP(17) + 0.2e1 * (-t163 * t206 + t194 * t239) * MDP(16) + (-t138 * t206 + t181 * t169) * t247 + (t139 * t206 - t181 * t170) * t246 + 0.2e1 * (-t109 * t206 + t133 * t124) * MDP(32) + (t110 * t206 + t133 * t125) * t243 + (-t118 * t143 + t119 * t142) * t245; (-t179 * t169 + t170 * t178) * MDP(19) + (t191 * t169 + t181 * t178) * MDP(23) + (-t191 * t170 + t181 * t179) * MDP(24) + (-t118 * t154 + t119 * t153 - t130 * t143 + t131 * t142) * MDP(25) + (t118 * t130 + t119 * t131 + t156 * t165) * MDP(26) + (-t137 * t124 - t125 * t136) * MDP(28) + (t140 * t124 + t133 * t136) * MDP(32) + (t140 * t125 + t133 * t137) * MDP(33) - t179 * t224 + t137 * t230 + (MDP(6) + (-t193 + t195) * MDP(12) + (-pkin(2) * t201 - t238) * MDP(16) + (-pkin(2) * t205 + t239) * MDP(17) - pkin(7) * MDP(9) + MDP(11) * t234) * t202 + (-pkin(7) * MDP(10) + MDP(7) - t254) * t206; MDP(8) + t193 * MDP(11) + 0.2e1 * t219 + t191 * t178 * t247 + (-t130 * t154 + t131 * t153) * t245 + (t130 ^ 2 + t131 ^ 2 + t165 ^ 2) * MDP(26) + 0.2e1 * t140 * t229 + 0.2e1 * (t205 * MDP(16) - t201 * MDP(17)) * pkin(2) + (MDP(18) * t179 + t178 * t248 + t191 * t246) * t179 + (MDP(27) * t137 + t136 * t244 + t140 * t243) * t137; t163 * MDP(16) - t164 * MDP(17) + (t173 * t142 - t171 * t143) * MDP(25) + (t118 * t171 + t119 * t173) * MDP(26) + t216 * t202 + (-MDP(15) - MDP(22) - t212 - t210) * t206 + t252; (t173 * t153 - t171 * t154) * MDP(25) + (t130 * t171 + t131 * t173) * MDP(26) + t254; (t171 ^ 2 + t173 ^ 2) * MDP(26) + 0.2e1 * t210 + 0.2e1 * t226 - 0.2e1 * t225 + t218; (-MDP(22) - t211) * t206 + ((t142 * t197 - t143 * t198) * MDP(25) + (t118 * t198 + t119 * t197) * MDP(26)) * pkin(4) + t252; ((t153 * t197 - t154 * t198) * MDP(25) + (t130 * t198 + t131 * t197) * MDP(26)) * pkin(4) + t209; (t147 + t184) * MDP(32) + (-t148 - t235) * MDP(33) + (t171 * t198 * MDP(26) + (t173 * MDP(26) - t199 * MDP(32) - t203 * MDP(33)) * t197) * pkin(4) + t210 + t221; (t197 ^ 2 + t198 ^ 2) * MDP(26) * pkin(4) ^ 2 + 0.2e1 * t223 - 0.2e1 * t222 + t221; t156 * MDP(26) + t124 * MDP(32) + t125 * MDP(33); t165 * MDP(26) + t137 * MDP(33) + t229; 0; 0; MDP(26); -t206 * MDP(31) + t251; t217; t212; t211; 0; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
