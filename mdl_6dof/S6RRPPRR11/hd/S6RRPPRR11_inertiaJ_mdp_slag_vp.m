% Calculate joint inertia matrix for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR11_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:42:47
% EndTime: 2019-03-09 09:42:51
% DurationCPUTime: 1.16s
% Computational Cost: add. (1235->210), mult. (2654->298), div. (0->0), fcn. (2877->10), ass. (0->96)
t174 = sin(pkin(6));
t182 = cos(qJ(2));
t217 = t174 * t182;
t180 = sin(qJ(2));
t218 = t174 * t180;
t157 = pkin(8) * t218;
t176 = cos(pkin(6));
t222 = pkin(1) * t182;
t203 = -pkin(2) - t222;
t131 = pkin(3) * t218 + t157 + (-qJ(4) + t203) * t176;
t177 = -pkin(2) - qJ(4);
t200 = -qJ(3) * t180 - pkin(1);
t135 = (t177 * t182 + t200) * t174;
t173 = sin(pkin(11));
t175 = cos(pkin(11));
t121 = t175 * t131 - t135 * t173;
t142 = -t173 * t217 + t175 * t176;
t117 = pkin(4) * t218 - pkin(9) * t142 + t121;
t122 = t173 * t131 + t175 * t135;
t141 = t176 * t173 + t175 * t217;
t118 = -pkin(9) * t141 + t122;
t179 = sin(qJ(5));
t223 = cos(qJ(5));
t114 = t223 * t117 - t179 * t118;
t115 = t179 * t117 + t223 * t118;
t127 = -t179 * t141 + t223 * t142;
t232 = t127 * MDP(21) + t114 * MDP(24) - t115 * MDP(25);
t178 = sin(qJ(6));
t181 = cos(qJ(6));
t190 = t178 * MDP(31) + t181 * MDP(32);
t231 = 0.2e1 * t176;
t149 = t173 * t179 - t223 * t175;
t229 = t149 * MDP(25);
t228 = -0.2e1 * pkin(2);
t227 = t149 ^ 2;
t150 = t223 * t173 + t179 * t175;
t146 = t150 ^ 2;
t225 = 0.2e1 * MDP(31);
t224 = 0.2e1 * MDP(32);
t221 = -pkin(9) + t177;
t126 = t223 * t141 + t142 * t179;
t220 = t126 * t150;
t152 = t221 * t173;
t201 = t221 * t175;
t133 = t152 * t179 - t223 * t201;
t219 = t133 * t149;
t216 = t178 * t181;
t162 = t173 * pkin(4) + qJ(3);
t144 = t176 * t180 * pkin(1) + pkin(8) * t217;
t156 = t173 ^ 2 + t175 ^ 2;
t215 = MDP(25) * t150;
t123 = t127 * t178 - t181 * t218;
t214 = t123 * MDP(29);
t124 = t127 * t181 + t178 * t218;
t213 = t124 * MDP(26);
t212 = t124 * MDP(28);
t211 = t126 * MDP(30);
t210 = t127 * MDP(20);
t208 = t127 * MDP(25);
t134 = t223 * t152 + t179 * t201;
t207 = t134 * MDP(25);
t164 = t176 * qJ(3);
t137 = -t164 - t144;
t136 = pkin(3) * t217 - t137;
t206 = t136 * MDP(18);
t202 = MDP(27) * t216;
t199 = t121 * t175 + t122 * t173;
t198 = -t144 * MDP(10) + (t176 * t222 - t157) * MDP(9);
t197 = t141 * MDP(15) + t142 * MDP(16);
t196 = MDP(15) * t175 - MDP(16) * t173;
t195 = t173 * MDP(15) + t175 * MDP(16);
t193 = t126 * MDP(24) + t208;
t192 = -MDP(28) * t181 + MDP(29) * t178;
t191 = MDP(31) * t181 - MDP(32) * t178;
t125 = pkin(4) * t141 + t136;
t189 = MDP(20) + t192;
t188 = MDP(24) + t191;
t187 = t178 * MDP(28) + t181 * MDP(29) - t190 * pkin(10);
t113 = pkin(10) * t218 + t115;
t116 = pkin(5) * t126 - pkin(10) * t127 + t125;
t110 = -t113 * t178 + t116 * t181;
t111 = t113 * t181 + t116 * t178;
t186 = t110 * MDP(31) - t111 * MDP(32) + t211 + t212 - t214;
t185 = -MDP(22) + t187;
t184 = t199 * MDP(18) + (-t141 * t173 - t142 * t175) * MDP(17);
t183 = qJ(3) ^ 2;
t172 = t181 ^ 2;
t170 = t178 ^ 2;
t148 = t156 * t177;
t140 = t203 * t176 + t157;
t138 = (-pkin(2) * t182 + t200) * t174;
t128 = pkin(5) * t150 + pkin(10) * t149 + t162;
t120 = t128 * t178 + t134 * t181;
t119 = t128 * t181 - t134 * t178;
t112 = -pkin(5) * t218 - t114;
t1 = [t127 ^ 2 * MDP(19) + (t121 ^ 2 + t122 ^ 2 + t136 ^ 2) * MDP(18) + (t137 ^ 2 + t138 ^ 2 + t140 ^ 2) * MDP(14) + t176 ^ 2 * MDP(8) + MDP(1) + (-0.2e1 * t123 * MDP(27) + t213) * t124 + (-0.2e1 * MDP(22) * t218 - 0.2e1 * t210 + t211 + 0.2e1 * t212 - 0.2e1 * t214) * t126 + 0.2e1 * (-t121 * t142 - t122 * t141) * MDP(17) + (t110 * t126 + t112 * t123) * t225 + (-t111 * t126 + t112 * t124) * t224 + (t140 * MDP(12) - t137 * MDP(13) + t198) * t231 + 0.2e1 * t197 * t136 + 0.2e1 * t193 * t125 + ((t180 * MDP(6) + t182 * MDP(7)) * t231 + 0.2e1 * (-t137 * MDP(11) + t138 * MDP(12)) * t182 + ((MDP(23) + MDP(4)) * t180 ^ 2 + 0.2e1 * (-MDP(10) * t180 + MDP(9) * t182) * pkin(1)) * t174 + 0.2e1 * (t140 * MDP(11) - t138 * MDP(13) + t121 * MDP(15) - t122 * MDP(16) + MDP(5) * t217 + t232) * t180) * t174; (-pkin(2) * t140 - qJ(3) * t137) * MDP(14) + (0.2e1 * t164 + t144) * MDP(13) + qJ(3) * t206 - t199 * MDP(17) + (t119 * t126 + t123 * t133) * MDP(31) + (qJ(3) * t141 + t136 * t173) * MDP(15) + (qJ(3) * t142 + t136 * t175) * MDP(16) + (-t120 * t126 + t124 * t133) * MDP(32) + t157 * MDP(12) + t184 * t177 + (MDP(8) + (t228 - t222) * MDP(12)) * t176 + t193 * t162 + (t125 * MDP(24) + t186 - t210) * t150 + (-t127 * MDP(19) - t181 * t213 - t125 * MDP(25) + (t123 * t181 + t124 * t178) * MDP(27) - t190 * t112 + t189 * t126) * t149 + ((qJ(3) * MDP(11) + MDP(7)) * t182 + (-pkin(2) * MDP(11) - t149 * MDP(21) - t150 * MDP(22) - t133 * MDP(24) + t196 * t177 + MDP(6) - t207) * t180) * t174 + t198; MDP(8) + MDP(12) * t228 + (pkin(2) ^ 2 + t183) * MDP(14) + (t156 * t177 ^ 2 + t183) * MDP(18) + t146 * MDP(30) + (t172 * MDP(26) + MDP(19) - 0.2e1 * t202) * t227 - 0.2e1 * t148 * MDP(17) + (t119 * t150 - t178 * t219) * t225 + (-t120 * t150 - t181 * t219) * t224 + 0.2e1 * (t150 * MDP(24) - t229) * t162 + 0.2e1 * (MDP(13) + t195) * qJ(3) + 0.2e1 * t150 * t189 * t149; t176 * MDP(12) + t140 * MDP(14) + (t123 * t149 - t178 * t220) * MDP(31) + (t124 * t149 - t181 * t220) * MDP(32) + (-MDP(24) * t149 + MDP(11) + t196 - t215) * t218 + t184; -pkin(2) * MDP(14) - t156 * MDP(17) + t148 * MDP(18) + MDP(12) + t190 * (-t146 - t227); MDP(18) * t156 + MDP(14); t188 * t126 + t197 + t206 + t208; qJ(3) * MDP(18) + t188 * t150 + t195 - t229; 0; MDP(18); MDP(23) * t218 + t178 * t213 + (-t123 * t178 + t124 * t181) * MDP(27) + (-pkin(5) * t123 - t112 * t181) * MDP(31) + (-pkin(5) * t124 + t112 * t178) * MDP(32) + t185 * t126 + t232; -t207 - t188 * t133 + t185 * t150 + (-MDP(21) - MDP(26) * t216 + (t170 - t172) * MDP(27) + t190 * pkin(5)) * t149; -t188 * t149 - t215; 0; MDP(26) * t170 + 0.2e1 * pkin(5) * t191 + MDP(23) + 0.2e1 * t202; t186; t150 * MDP(30) + t119 * MDP(31) - t120 * MDP(32) + t192 * t149; -t190 * t150; t191; t187; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
