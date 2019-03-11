% Calculate joint inertia matrix for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR3_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:51:35
% EndTime: 2019-03-09 00:51:37
% DurationCPUTime: 0.74s
% Computational Cost: add. (914->180), mult. (1953->260), div. (0->0), fcn. (2181->12), ass. (0->97)
t171 = sin(qJ(4));
t176 = cos(qJ(4));
t219 = pkin(9) + pkin(10);
t155 = t219 * t171;
t156 = t219 * t176;
t170 = sin(qJ(5));
t175 = cos(qJ(5));
t132 = -t175 * t155 - t170 * t156;
t133 = -t170 * t155 + t175 * t156;
t150 = t170 * t171 - t175 * t176;
t151 = t170 * t176 + t175 * t171;
t118 = -t151 * pkin(11) + t132;
t119 = -t150 * pkin(11) + t133;
t169 = sin(qJ(6));
t174 = cos(qJ(6));
t123 = t174 * t150 + t169 * t151;
t124 = -t169 * t150 + t174 * t151;
t189 = t124 * MDP(28) - t123 * MDP(29) + (t174 * t118 - t169 * t119) * MDP(31) - (t169 * t118 + t174 * t119) * MDP(32);
t180 = t151 * MDP(21) - t150 * MDP(22) + t132 * MDP(24) - t133 * MDP(25) + t189;
t230 = t180 - (t171 * MDP(17) + t176 * MDP(18)) * pkin(9) + t171 * MDP(14) + t176 * MDP(15);
t196 = t169 * MDP(32);
t181 = (t174 * MDP(31) - t196) * pkin(5);
t195 = t175 * MDP(24);
t229 = (-t170 * MDP(25) + t195) * pkin(4);
t172 = sin(qJ(3));
t140 = t151 * t172;
t141 = t150 * t172;
t116 = t174 * t140 - t169 * t141;
t117 = -t169 * t140 - t174 * t141;
t206 = t117 * MDP(28) - t116 * MDP(29);
t228 = -t141 * MDP(21) - t140 * MDP(22);
t177 = cos(qJ(3));
t154 = -t177 * pkin(3) - t172 * pkin(9) - pkin(2);
t149 = t176 * t154;
t209 = t172 * t176;
t217 = pkin(8) * t171;
t125 = -pkin(10) * t209 + t149 + (-pkin(4) - t217) * t177;
t215 = pkin(8) * t177;
t193 = t176 * t215;
t129 = t193 + (-pkin(10) * t172 + t154) * t171;
t110 = t175 * t125 - t170 * t129;
t111 = t170 * t125 + t175 * t129;
t102 = -t177 * pkin(5) + t141 * pkin(11) + t110;
t107 = -t140 * pkin(11) + t111;
t95 = t174 * t102 - t169 * t107;
t96 = t169 * t102 + t174 * t107;
t226 = t95 * MDP(31) - t96 * MDP(32) + t206;
t227 = t110 * MDP(24) - t111 * MDP(25) + t226 + t228;
t223 = -2 * MDP(20);
t222 = 0.2e1 * MDP(25);
t221 = -2 * MDP(27);
t220 = 0.2e1 * MDP(32);
t218 = pkin(4) * t170;
t216 = pkin(8) * t176;
t168 = cos(pkin(6));
t167 = sin(pkin(6));
t173 = sin(qJ(2));
t213 = t167 * t173;
t145 = t168 * t172 + t177 * t213;
t178 = cos(qJ(2));
t212 = t167 * t178;
t127 = -t145 * t171 - t176 * t212;
t128 = t145 * t176 - t171 * t212;
t112 = t175 * t127 - t170 * t128;
t113 = t170 * t127 + t175 * t128;
t100 = t169 * t112 + t174 * t113;
t99 = t174 * t112 - t169 * t113;
t214 = t99 * MDP(31) - t100 * MDP(32);
t211 = t171 * t172;
t210 = t171 * t176;
t153 = pkin(4) * t211 + t172 * pkin(8);
t205 = MDP(11) * t172;
t202 = t123 * MDP(31);
t201 = t124 * MDP(26);
t160 = t175 * pkin(4) + pkin(5);
t157 = t174 * t160;
t200 = (-t169 * t218 + t157) * MDP(31);
t199 = (t169 * t160 + t174 * t218) * MDP(32);
t198 = t150 * MDP(24);
t197 = t151 * MDP(19);
t194 = MDP(23) + MDP(30);
t161 = -t176 * pkin(4) - pkin(3);
t192 = MDP(13) * t210;
t191 = MDP(16) + t194;
t190 = t112 * MDP(24) - t113 * MDP(25) + t214;
t188 = t176 * MDP(14) - t171 * MDP(15);
t186 = t176 * MDP(17) - t171 * MDP(18);
t182 = MDP(30) - t199 + t200;
t165 = t176 ^ 2;
t164 = t172 ^ 2;
t163 = t171 ^ 2;
t144 = -t168 * t177 + t172 * t213;
t137 = t150 * pkin(5) + t161;
t136 = t171 * t154 + t193;
t135 = -t171 * t215 + t149;
t126 = t140 * pkin(5) + t153;
t1 = [MDP(1); (-t127 * t177 + t144 * t211) * MDP(17) + (t128 * t177 + t144 * t209) * MDP(18) + (-t112 * t177 + t144 * t140) * MDP(24) + (t113 * t177 - t144 * t141) * MDP(25) + (t144 * t116 - t99 * t177) * MDP(31) + (t100 * t177 + t144 * t117) * MDP(32) + (-t173 * MDP(4) + (MDP(10) * t177 + MDP(3) - t205) * t178) * t167; -0.2e1 * pkin(2) * t205 + MDP(2) - (-t141 * MDP(19) + t140 * t223) * t141 + (t117 * MDP(26) + t116 * t221) * t117 + t191 * t177 ^ 2 + (t165 * MDP(12) + MDP(5) - 0.2e1 * t192) * t164 + 0.2e1 * (pkin(2) * MDP(10) + (MDP(6) - t188) * t172 - t228 - t206) * t177 + 0.2e1 * (t126 * t116 - t95 * t177) * MDP(31) + (t126 * t117 + t96 * t177) * t220 + 0.2e1 * (-t110 * t177 + t153 * t140) * MDP(24) + (t111 * t177 - t153 * t141) * t222 + 0.2e1 * (-t135 * t177 + t164 * t217) * MDP(17) + 0.2e1 * (t136 * t177 + t164 * t216) * MDP(18); -t145 * MDP(11) + (t151 * MDP(25) + t124 * MDP(32) - MDP(10) - t186 + t198 + t202) * t144; -t141 * t197 + (-t151 * t140 + t141 * t150) * MDP(20) + (t161 * t140 + t153 * t150) * MDP(24) + (-t161 * t141 + t153 * t151) * MDP(25) + t117 * t201 + (-t124 * t116 - t117 * t123) * MDP(27) + (t137 * t116 + t126 * t123) * MDP(31) + (t137 * t117 + t126 * t124) * MDP(32) + (MDP(7) - pkin(8) * MDP(10) + MDP(12) * t210 + (-t163 + t165) * MDP(13) + (-pkin(3) * t171 - t216) * MDP(17) + (-pkin(3) * t176 + t217) * MDP(18)) * t172 + (-pkin(8) * MDP(11) + MDP(8) - t230) * t177; 0.2e1 * t192 + 0.2e1 * t161 * t198 + 0.2e1 * t137 * t202 + t163 * MDP(12) + MDP(9) + 0.2e1 * t186 * pkin(3) + (t150 * t223 + t161 * t222 + t197) * t151 + (t123 * t221 + t137 * t220 + t201) * t124; t127 * MDP(17) - t128 * MDP(18) + t190; t135 * MDP(17) - t136 * MDP(18) + t188 * t172 + (-MDP(16) - MDP(23) - t182 - t229) * t177 + t227; t230; t191 - 0.2e1 * t199 + 0.2e1 * t200 + 0.2e1 * t229; t190; (-t194 - t181) * t177 + t227; t180; (t174 * pkin(5) + t157) * MDP(31) + (-pkin(5) - t160) * t196 + (t195 + (-MDP(31) * t169 - MDP(32) * t174 - MDP(25)) * t170) * pkin(4) + t194; 0.2e1 * t181 + t194; t214; -t177 * MDP(30) + t226; t189; t182; MDP(30) + t181; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
