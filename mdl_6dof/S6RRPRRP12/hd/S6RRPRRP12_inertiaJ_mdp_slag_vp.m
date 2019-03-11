% Calculate joint inertia matrix for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRRP12_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:53:35
% EndTime: 2019-03-09 12:53:38
% DurationCPUTime: 1.05s
% Computational Cost: add. (960->197), mult. (1583->265), div. (0->0), fcn. (1513->6), ass. (0->81)
t154 = sin(qJ(4));
t156 = cos(qJ(4));
t158 = -pkin(2) - pkin(8);
t197 = -pkin(9) + t158;
t132 = t197 * t154;
t153 = sin(qJ(5));
t177 = t197 * t156;
t199 = cos(qJ(5));
t112 = t153 * t132 - t199 * t177;
t113 = t199 * t132 + t153 * t177;
t180 = t199 * t156;
t128 = t153 * t154 - t180;
t164 = -t153 * t156 - t199 * t154;
t183 = -MDP(27) - MDP(29);
t207 = MDP(28) - MDP(31);
t165 = -t128 * MDP(24) + MDP(25) * t164 + t183 * t112 - t207 * t113;
t212 = (t158 * MDP(20) + MDP(17)) * t156 + (-t158 * MDP(21) - MDP(18)) * t154 + t165;
t171 = t183 * t128 + t207 * t164;
t210 = pkin(3) + pkin(7);
t206 = -2 * MDP(23);
t157 = cos(qJ(2));
t193 = t154 * t157;
t116 = t153 * t193 - t157 * t180;
t117 = t164 * t157;
t205 = t117 * MDP(24) + t116 * MDP(25);
t155 = sin(qJ(2));
t134 = t210 * t155;
t131 = t156 * t134;
t176 = -t155 * qJ(3) - pkin(1);
t126 = t158 * t157 + t176;
t175 = pkin(9) * t157 - t126;
t102 = t155 * pkin(4) + t175 * t154 + t131;
t195 = t154 * t134;
t104 = -t175 * t156 + t195;
t98 = t153 * t102 + t199 * t104;
t203 = t155 * MDP(26) - t98 * MDP(28) + t205;
t202 = 0.2e1 * t157;
t201 = 2 * MDP(29);
t200 = 2 * MDP(31);
t198 = t155 * pkin(5);
t196 = pkin(2) * MDP(14);
t194 = t154 * t156;
t192 = t156 * t157;
t139 = t154 * pkin(4) + qJ(3);
t135 = t210 * t157;
t190 = t117 * MDP(22);
t145 = t153 * pkin(4);
t138 = t145 + qJ(6);
t188 = t138 * MDP(31);
t187 = t154 * MDP(20);
t186 = t156 * MDP(21);
t185 = pkin(7) ^ 2 * MDP(14);
t184 = MDP(19) + MDP(26);
t144 = t155 * qJ(6);
t95 = t144 + t98;
t119 = pkin(4) * t192 + t135;
t179 = MDP(16) * t194;
t178 = t128 ^ 2 + t164 ^ 2;
t174 = MDP(12) - t196;
t173 = pkin(5) * t201 + MDP(26);
t172 = -t199 * t102 + t153 * t104;
t170 = t128 * pkin(5) + qJ(6) * t164;
t96 = t172 - t198;
t169 = t96 * t128 - t164 * t95;
t141 = t199 * pkin(4) + pkin(5);
t168 = t128 * t141 + t138 * t164;
t167 = -t154 * MDP(17) - t156 * MDP(18);
t166 = t156 * MDP(20) - t154 * MDP(21);
t163 = t199 * MDP(27) - t153 * MDP(28);
t162 = t163 * pkin(4);
t161 = pkin(7) * MDP(14) + MDP(11) + t166;
t152 = t157 ^ 2;
t151 = t156 ^ 2;
t150 = t155 ^ 2;
t149 = t154 ^ 2;
t133 = -t157 * pkin(2) + t176;
t107 = t156 * t126 + t195;
t106 = -t154 * t126 + t131;
t105 = -pkin(5) * t164 + t128 * qJ(6) + t139;
t99 = -t116 * pkin(5) - t117 * qJ(6) + t119;
t1 = [MDP(1) + (t95 ^ 2 + t96 ^ 2 + t99 ^ 2) * MDP(32) + pkin(1) * MDP(9) * t202 + (MDP(12) * t202 + MDP(14) * t133) * t133 + (t149 * MDP(15) + 0.2e1 * t179 + t185) * t152 - (t116 * t206 - t190) * t117 + (MDP(4) + t184 + t185) * t150 + 0.2e1 * (-pkin(1) * MDP(10) - t133 * MDP(13) + (MDP(5) + t167) * t157 + t205) * t155 + 0.2e1 * (t95 * t116 + t96 * t117) * MDP(30) + (-t99 * t116 - t96 * t155) * t201 + (-t99 * t117 + t95 * t155) * t200 + 0.2e1 * (-t119 * t116 - t155 * t172) * MDP(27) + 0.2e1 * (t119 * t117 - t98 * t155) * MDP(28) + 0.2e1 * (-t107 * t155 - t135 * t193) * MDP(21) + 0.2e1 * (t106 * t155 + t135 * t192) * MDP(20) + 0.2e1 * (t150 + t152) * MDP(11) * pkin(7); (-t128 * t116 + t117 * t164) * MDP(23) + (-t139 * t116 - t119 * t164) * MDP(27) + (t139 * t117 - t119 * t128) * MDP(28) + (-t105 * t116 - t164 * t99) * MDP(29) + (t112 * t117 + t113 * t116 - t169) * MDP(30) + (-t105 * t117 + t99 * t128) * MDP(31) + (t99 * t105 + t96 * t112 + t95 * t113) * MDP(32) - t128 * t190 + (t186 + t187) * t135 + (MDP(7) + (t149 - t151) * MDP(16) - MDP(15) * t194 + (MDP(13) - MDP(10)) * pkin(7) + t161 * qJ(3)) * t157 + (-pkin(2) * MDP(11) + MDP(6) + (-MDP(9) + t174) * pkin(7) + t212) * t155; MDP(8) + t151 * MDP(15) - 0.2e1 * t179 + (t105 ^ 2 + t112 ^ 2 + t113 ^ 2) * MDP(32) + (-0.2e1 * MDP(12) + t196) * pkin(2) - 0.2e1 * (t139 * MDP(27) + t105 * MDP(29) - t113 * MDP(30)) * t164 + (MDP(22) * t128 - 0.2e1 * t139 * MDP(28) - 0.2e1 * t112 * MDP(30) + t105 * t200 + t164 * t206) * t128 + (MDP(14) * qJ(3) + 0.2e1 * MDP(13) + 0.2e1 * t186 + 0.2e1 * t187) * qJ(3); (-t116 * t164 + t128 * t117) * MDP(30) + t169 * MDP(32) + (t161 + t171) * t155; -t178 * MDP(30) + (t112 * t128 - t113 * t164) * MDP(32) + t174; t178 * MDP(32) + MDP(14); t106 * MDP(20) - t107 * MDP(21) + (t138 * t116 - t141 * t117) * MDP(30) + t95 * MDP(31) + (t95 * t138 - t96 * t141) * MDP(32) + t167 * t157 + (MDP(19) + (pkin(5) + t141) * MDP(29) + t188 + t162) * t155 + t183 * t172 + t203; t168 * MDP(30) + (-t112 * t141 + t113 * t138) * MDP(32) + t212; -t168 * MDP(32) + t166 + t171; (t138 ^ 2 + t141 ^ 2) * MDP(32) + 0.2e1 * t162 + t141 * t201 + 0.2e1 * t188 + t184; -t172 * MDP(27) + (-t172 + 0.2e1 * t198) * MDP(29) + (-pkin(5) * t117 + t116 * qJ(6)) * MDP(30) + (0.2e1 * t144 + t98) * MDP(31) + (-t96 * pkin(5) + t95 * qJ(6)) * MDP(32) + t203; t170 * MDP(30) + (-t112 * pkin(5) + t113 * qJ(6)) * MDP(32) + t165; -t170 * MDP(32) + t171; (0.2e1 * qJ(6) + t145) * MDP(31) + (t141 * pkin(5) + t138 * qJ(6)) * MDP(32) + (t199 * MDP(29) + t163) * pkin(4) + t173; qJ(6) * t200 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(32) + t173; -t155 * MDP(29) + t117 * MDP(30) + t96 * MDP(32); -t128 * MDP(30) + t112 * MDP(32); t128 * MDP(32); -t141 * MDP(32) - MDP(29); -MDP(32) * pkin(5) - MDP(29); MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
