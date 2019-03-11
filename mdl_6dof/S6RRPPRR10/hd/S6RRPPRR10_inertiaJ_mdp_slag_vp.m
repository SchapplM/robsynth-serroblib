% Calculate joint inertia matrix for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR10_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:37:03
% EndTime: 2019-03-09 09:37:05
% DurationCPUTime: 0.58s
% Computational Cost: add. (816->157), mult. (1479->210), div. (0->0), fcn. (1538->8), ass. (0->81)
t214 = pkin(3) + pkin(7);
t166 = -pkin(2) - qJ(4);
t172 = cos(qJ(2));
t169 = sin(qJ(2));
t191 = -qJ(3) * t169 - pkin(1);
t140 = t166 * t172 + t191;
t149 = t214 * t169;
t165 = cos(pkin(10));
t145 = t165 * t149;
t164 = sin(pkin(10));
t116 = pkin(4) * t169 + t145 + (pkin(8) * t172 - t140) * t164;
t128 = t165 * t140 + t164 * t149;
t200 = t165 * t172;
t117 = -pkin(8) * t200 + t128;
t168 = sin(qJ(5));
t171 = cos(qJ(5));
t106 = t171 * t116 - t117 * t168;
t107 = t116 * t168 + t117 * t171;
t201 = t164 * t168;
t133 = -t171 * t200 + t172 * t201;
t142 = t164 * t171 + t165 * t168;
t134 = t142 * t172;
t213 = -t134 * MDP(21) + t133 * MDP(22) + t106 * MDP(24) - t107 * MDP(25);
t167 = sin(qJ(6));
t170 = cos(qJ(6));
t114 = -t170 * t133 - t134 * t167;
t115 = t133 * t167 - t134 * t170;
t210 = t115 * MDP(28) - t114 * MDP(29);
t203 = -pkin(8) + t166;
t146 = t203 * t164;
t147 = t203 * t165;
t129 = -t146 * t168 + t171 * t147;
t130 = t146 * t171 + t147 * t168;
t143 = t165 * t171 - t201;
t112 = -pkin(9) * t143 + t129;
t113 = -pkin(9) * t142 + t130;
t186 = t170 * t142 + t143 * t167;
t189 = -t142 * t167 + t170 * t143;
t188 = t189 * MDP(28) - t186 * MDP(29) + (t112 * t170 - t113 * t167) * MDP(31) - (t112 * t167 + t113 * t170) * MDP(32);
t209 = t143 * MDP(21) - t142 * MDP(22) + t129 * MDP(24) - t130 * MDP(25) + t188;
t199 = MDP(31) * t189 - MDP(32) * t186;
t208 = MDP(24) * t143 - MDP(25) * t142 + t199;
t154 = t164 * pkin(4) + qJ(3);
t132 = pkin(5) * t142 + t154;
t207 = 0.2e1 * t132;
t206 = 0.2e1 * t154;
t205 = -2 * MDP(27);
t204 = pkin(5) * t169;
t101 = pkin(9) * t133 + t107;
t202 = t101 * t170;
t150 = t214 * t172;
t152 = t164 ^ 2 + t165 ^ 2;
t127 = -t140 * t164 + t145;
t108 = t127 * t165 + t128 * t164;
t198 = t108 * MDP(18);
t197 = t115 * MDP(26);
t196 = t186 * MDP(31);
t195 = t134 * MDP(19);
t194 = t142 * MDP(24);
t193 = MDP(23) + MDP(30);
t192 = t169 * MDP(30) + t210;
t138 = pkin(4) * t200 + t150;
t100 = pkin(9) * t134 + t106 + t204;
t97 = t170 * t100 - t101 * t167;
t190 = -pkin(2) * MDP(14) + MDP(12);
t187 = t97 * MDP(31) - (t100 * t167 + t202) * MDP(32);
t185 = MDP(15) * t165 - MDP(16) * t164;
t184 = MDP(15) * t164 + MDP(16) * t165;
t181 = -t133 * MDP(24) - t134 * MDP(25);
t179 = t114 * MDP(31) + t115 * MDP(32);
t178 = (MDP(31) * t170 - MDP(32) * t167) * pkin(5);
t177 = qJ(3) * MDP(18) + t184;
t176 = MDP(14) * pkin(7) + MDP(11) + t185;
t174 = pkin(7) ^ 2;
t173 = qJ(3) ^ 2;
t163 = t172 ^ 2;
t162 = t169 ^ 2;
t148 = -pkin(2) * t172 + t191;
t141 = t152 * t166;
t118 = -pkin(5) * t133 + t138;
t1 = [(t127 ^ 2 + t128 ^ 2 + t150 ^ 2) * MDP(18) + (t148 ^ 2 + t163 * t174) * MDP(14) + MDP(1) - (0.2e1 * t133 * MDP(20) - t195) * t134 + (t114 * t205 + t197) * t115 + (t174 * MDP(14) + MDP(4) + t193) * t162 + 0.2e1 * t181 * t138 + 0.2e1 * t179 * t118 + 0.2e1 * (t162 + t163) * MDP(11) * pkin(7) + 0.2e1 * (-pkin(1) * MDP(10) - t148 * MDP(13) + t127 * MDP(15) - t128 * MDP(16) + t172 * MDP(5) + t187 + t210 + t213) * t169 + 0.2e1 * ((t127 * t164 - t128 * t165) * MDP(17) + t185 * t150 + t148 * MDP(12) + pkin(1) * MDP(9)) * t172; -t108 * MDP(17) + (t133 * t143 + t134 * t142) * MDP(20) + (-t133 * t154 + t138 * t142) * MDP(24) + (-t134 * t154 + t138 * t143) * MDP(25) + (-t114 * t189 - t115 * t186) * MDP(27) + (t114 * t132 + t118 * t186) * MDP(31) + (t115 * t132 + t118 * t189) * MDP(32) - t143 * t195 + t189 * t197 + t166 * t198 + t177 * t150 + (MDP(7) + (-MDP(10) + MDP(13)) * pkin(7) + t176 * qJ(3)) * t172 + (-pkin(2) * MDP(11) + MDP(6) + t185 * t166 + (-MDP(9) + t190) * pkin(7) + t209) * t169; MDP(8) - 0.2e1 * pkin(2) * MDP(12) + (pkin(2) ^ 2 + t173) * MDP(14) - 0.2e1 * t141 * MDP(17) + (t152 * t166 ^ 2 + t173) * MDP(18) + t194 * t206 + t196 * t207 + (MDP(19) * t143 - 0.2e1 * t142 * MDP(20) + MDP(25) * t206) * t143 + (MDP(26) * t189 + MDP(32) * t207 + t186 * t205) * t189 + 0.2e1 * (MDP(13) + t184) * qJ(3); t198 + (t176 + t208) * t169; -MDP(17) * t152 + MDP(18) * t141 + t190; MDP(18) * t152 + MDP(14); t150 * MDP(18) + t185 * t172 + t179 + t181; t143 * MDP(25) + MDP(32) * t189 + t177 + t194 + t196; 0; MDP(18); t169 * MDP(23) + (t170 * t204 + t97) * MDP(31) + (-t202 + (-t100 - t204) * t167) * MDP(32) + t192 + t213; t209; t208; 0; 0.2e1 * t178 + t193; t187 + t192; t188; t199; 0; MDP(30) + t178; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
