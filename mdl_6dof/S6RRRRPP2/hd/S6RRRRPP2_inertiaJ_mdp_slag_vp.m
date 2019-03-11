% Calculate joint inertia matrix for
% S6RRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP2_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:51:42
% EndTime: 2019-03-09 20:51:46
% DurationCPUTime: 1.10s
% Computational Cost: add. (1252->247), mult. (2152->302), div. (0->0), fcn. (2181->6), ass. (0->92)
t182 = sin(qJ(3));
t164 = pkin(2) * t182 + pkin(9);
t181 = sin(qJ(4));
t179 = t181 ^ 2;
t184 = cos(qJ(4));
t180 = t184 ^ 2;
t216 = t179 + t180;
t219 = t216 * t164;
t241 = t181 * MDP(20) + t184 * MDP(21);
t183 = sin(qJ(2));
t230 = cos(qJ(3));
t231 = cos(qJ(2));
t150 = t182 * t183 - t230 * t231;
t240 = 0.2e1 * t150;
t206 = 0.2e1 * t181;
t205 = 0.2e1 * t184;
t239 = pkin(8) + pkin(7);
t238 = t184 * pkin(4) + t181 * qJ(5);
t151 = t182 * t231 + t230 * t183;
t166 = -t231 * pkin(2) - pkin(1);
t126 = t150 * pkin(3) - t151 * pkin(9) + t166;
t157 = t239 * t183;
t158 = t239 * t231;
t133 = -t182 * t157 + t230 * t158;
t120 = t184 * t126 - t181 * t133;
t116 = -t150 * pkin(4) - t120;
t121 = t181 * t126 + t184 * t133;
t237 = t120 * MDP(23) - t121 * MDP(24) - t116 * MDP(25);
t132 = t230 * t157 + t158 * t182;
t236 = t132 * MDP(24);
t185 = pkin(4) + pkin(5);
t225 = qJ(5) * t184;
t235 = t181 * t185 - t225;
t234 = 0.2e1 * t166;
t233 = 2 * MDP(26);
t232 = 2 * MDP(31);
t229 = pkin(9) * t150;
t165 = -t230 * pkin(2) - pkin(3);
t228 = pkin(3) - t165;
t227 = MDP(28) * pkin(9);
t226 = pkin(4) * MDP(25);
t224 = t150 * t164;
t223 = t151 * t184;
t173 = t181 * qJ(6);
t222 = t181 * t184;
t142 = t165 - t238;
t177 = t184 * pkin(5);
t135 = t177 - t142;
t204 = pkin(3) + t238;
t143 = t177 + t204;
t221 = t135 + t143;
t220 = -t142 + t204;
t218 = t184 * MDP(29) + t181 * MDP(30);
t217 = t216 * pkin(9);
t215 = MDP(28) * t116;
t197 = pkin(4) * t181 - t225;
t122 = t197 * t151 + t132;
t214 = MDP(28) * t122;
t213 = MDP(28) * t164;
t212 = MDP(29) * t151;
t211 = MDP(30) * t151;
t209 = t122 * MDP(27);
t208 = 0.2e1 * t231;
t207 = MDP(26) - MDP(31);
t144 = t150 * qJ(5);
t115 = t144 + t121;
t203 = 0.2e1 * t144 + t121;
t202 = MDP(19) * t222;
t201 = t179 * MDP(18) + MDP(15) + 0.2e1 * t202;
t200 = -pkin(4) * MDP(28) - MDP(25);
t199 = -t197 * MDP(26) + MDP(31) * t235 + t241;
t198 = -pkin(3) * t151 - t229;
t196 = t151 * t204 + t229;
t195 = t142 * t151 - t224;
t194 = t151 * t165 - t224;
t192 = -t132 * MDP(23) - t122 * MDP(25);
t191 = (t230 * MDP(16) - t182 * MDP(17)) * pkin(2);
t190 = t181 * t236 + (t115 * t184 + t116 * t181) * MDP(26) - t132 * MDP(16) - t133 * MDP(17) + ((-t179 + t180) * MDP(19) + MDP(18) * t222 + MDP(13)) * t151 + (-MDP(14) + t241) * t150;
t189 = (MDP(28) * qJ(5) - MDP(24) + MDP(27)) * t184 + (-MDP(23) + t200) * t181;
t187 = qJ(5) ^ 2;
t169 = t181 * MDP(26);
t156 = (pkin(9) - qJ(6)) * t184;
t155 = pkin(9) * t181 - t173;
t146 = (-qJ(6) + t164) * t184;
t145 = t164 * t181 - t173;
t138 = t151 * t173;
t119 = -t151 * t235 - t132;
t118 = t119 * t184;
t117 = t119 * t181;
t114 = t138 + t115;
t113 = -pkin(5) * t150 - qJ(6) * t223 + t116;
t1 = [pkin(1) * MDP(9) * t208 + (t115 ^ 2 + t116 ^ 2 + t122 ^ 2) * MDP(28) + (t113 ^ 2 + t114 ^ 2 + t119 ^ 2) * MDP(32) + MDP(1) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t183 + MDP(5) * t208) * t183 + (MDP(16) * t234 + t150 * MDP(22)) * t150 + (t115 * MDP(27) - t113 * MDP(29) + t114 * MDP(30) + t237) * t240 + (MDP(17) * t234 + (t184 * MDP(20) - t181 * MDP(21) - MDP(12)) * t240 + (t116 * MDP(26) + t119 * MDP(30) - t113 * MDP(31) - t209 + t236) * t205 + (-t115 * MDP(26) - t119 * MDP(29) + t114 * MDP(31) - t192) * t206 + (t180 * MDP(18) + MDP(11) - 0.2e1 * t202) * t151) * t151; (-t231 * MDP(10) - t183 * MDP(9)) * pkin(7) + (t194 * MDP(23) + t195 * MDP(25) - t209 + t116 * t213 - t135 * t212 + (t146 * t151 - t113) * MDP(31)) * t181 + (-t145 * t150 + t118) * MDP(29) + (t146 * t150 + t117) * MDP(30) + t183 * MDP(6) + t142 * t214 + t190 + (t113 * t145 + t114 * t146 + t119 * t135) * MDP(32) + t231 * MDP(7) + (t194 * MDP(24) - t195 * MDP(27) + t115 * t213 + t135 * t211 + (-t145 * t151 - t114) * MDP(31) + t192) * t184; MDP(8) + (t216 * t164 ^ 2 + t142 ^ 2) * MDP(28) + (t135 ^ 2 + t145 ^ 2 + t146 ^ 2) * MDP(32) + 0.2e1 * t191 + (-MDP(23) * t165 - MDP(25) * t142 + MDP(29) * t135) * t205 + (MDP(24) * t165 - MDP(27) * t142 + MDP(30) * t135) * t206 + t219 * t233 + (-t145 * t181 - t146 * t184) * t232 + t201; (t113 * t155 + t114 * t156 + t119 * t143) * MDP(32) + (-t150 * t155 + t118) * MDP(29) + (t150 * t156 + t117) * MDP(30) + t190 + (t198 * MDP(23) - t196 * MDP(25) - t209 + pkin(9) * t215 - t143 * t212 + (t151 * t156 - t113) * MDP(31)) * t181 + (t198 * MDP(24) + t196 * MDP(27) + t115 * t227 + t143 * t211 + (-t151 * t155 - t114) * MDP(31) + t192) * t184 - t204 * t214; (t217 + t219) * MDP(26) + (pkin(9) * t219 - t142 * t204) * MDP(28) + (t135 * t143 + t145 * t155 + t146 * t156) * MDP(32) + t191 + (t228 * MDP(23) + t220 * MDP(25) + t221 * MDP(29) + (-t146 - t156) * MDP(31)) * t184 + (-t228 * MDP(24) + t220 * MDP(27) + t221 * MDP(30) + (-t145 - t155) * MDP(31)) * t181 + t201; (t216 * pkin(9) ^ 2 + t204 ^ 2) * MDP(28) + (t143 ^ 2 + t155 ^ 2 + t156 ^ 2) * MDP(32) + (MDP(23) * pkin(3) + MDP(25) * t204 + MDP(29) * t143) * t205 + (-MDP(24) * pkin(3) + MDP(27) * t204 + MDP(30) * t143) * t206 + t217 * t233 + (-t155 * t181 - t156 * t184) * t232 + t201; t203 * MDP(27) + (-pkin(4) * t116 + qJ(5) * t115) * MDP(28) - t116 * MDP(29) + (t138 + t203) * MDP(30) + (qJ(5) * t114 - t113 * t185) * MDP(32) + (MDP(22) + t226 + (pkin(5) + t185) * MDP(29)) * t150 + ((-t207 * qJ(5) - MDP(21)) * t181 + (-MDP(26) * pkin(4) + MDP(29) * qJ(6) + MDP(31) * t185 + MDP(20)) * t184) * t151 + t237; -t145 * MDP(29) + t146 * MDP(30) + (qJ(5) * t146 - t145 * t185) * MDP(32) + t189 * t164 + t199; -t155 * MDP(29) + t156 * MDP(30) + (qJ(5) * t156 - t155 * t185) * MDP(32) + t189 * pkin(9) + t199; MDP(22) + 0.2e1 * t226 + (pkin(4) ^ 2 + t187) * MDP(28) + 0.2e1 * t185 * MDP(29) + (t185 ^ 2 + t187) * MDP(32) + 0.2e1 * (MDP(27) + MDP(30)) * qJ(5); t215 + MDP(32) * t113 + t207 * t223 + (-MDP(25) - MDP(29)) * t150; MDP(32) * t145 + t169 + (-MDP(31) + t213) * t181; MDP(32) * t155 + t169 + (-MDP(31) + t227) * t181; -MDP(32) * t185 - MDP(29) + t200; MDP(28) + MDP(32); MDP(32) * t119 + (-t181 * MDP(29) + t184 * MDP(30)) * t151; MDP(32) * t135 + t218; MDP(32) * t143 + t218; 0; 0; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
