% Calculate joint inertia matrix for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRRRP8_inertiaJ_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:24:35
% EndTime: 2019-03-09 06:24:38
% DurationCPUTime: 0.65s
% Computational Cost: add. (857->163), mult. (1403->215), div. (0->0), fcn. (1419->6), ass. (0->81)
t148 = sin(qJ(4));
t136 = t148 * pkin(3) + pkin(9);
t147 = sin(qJ(5));
t145 = t147 ^ 2;
t150 = cos(qJ(5));
t146 = t150 ^ 2;
t184 = t145 + t146;
t186 = t184 * t136;
t149 = sin(qJ(3));
t151 = cos(qJ(3));
t200 = cos(qJ(4));
t124 = t148 * t151 + t149 * t200;
t121 = t124 ^ 2;
t123 = t148 * t149 - t151 * t200;
t203 = t123 ^ 2;
t211 = t121 + t203;
t209 = t147 * MDP(23) + t150 * MDP(24);
t179 = t147 * MDP(27);
t208 = t179 - MDP(19);
t207 = MDP(26) + MDP(28);
t206 = -MDP(27) + MDP(30);
t201 = 2 * MDP(29);
t199 = (pkin(1) * MDP(6));
t198 = pkin(9) * t124;
t197 = t124 * pkin(5);
t177 = t200 * pkin(3);
t137 = -t177 - pkin(4);
t196 = pkin(4) - t137;
t152 = -pkin(1) - pkin(7);
t195 = -pkin(8) + t152;
t134 = t149 * pkin(3) + qJ(2);
t104 = t124 * pkin(4) + t123 * pkin(9) + t134;
t127 = t195 * t149;
t128 = t195 * t151;
t110 = t127 * t200 + t148 * t128;
t172 = -t150 * t104 + t147 * t110;
t96 = t172 - t197;
t194 = t96 * t147;
t193 = t124 * qJ(6);
t192 = t124 * t136;
t191 = t147 * t150;
t98 = t147 * t104 + t150 * t110;
t95 = t98 + t193;
t190 = t95 * MDP(31);
t189 = t96 * MDP(31);
t109 = t148 * t127 - t128 * t200;
t166 = -pkin(5) * t147 + t150 * qJ(6);
t99 = t123 * t166 + t109;
t188 = t99 * MDP(30);
t167 = t150 * pkin(5) + t147 * qJ(6);
t129 = -pkin(4) - t167;
t120 = -t177 + t129;
t187 = -t120 - t129;
t185 = pkin(9) * t184;
t183 = MDP(31) * t136;
t182 = MDP(31) * t147;
t181 = t120 * MDP(31);
t180 = t129 * MDP(31);
t178 = t147 * MDP(30);
t176 = t166 * MDP(29) + t209;
t175 = MDP(22) * t191;
t174 = t145 * MDP(21) + MDP(18) + 0.2e1 * t175;
t173 = -MDP(31) * pkin(5) - MDP(28);
t171 = t184 * MDP(31);
t170 = (MDP(29) * t184 - MDP(20)) * t124 + t208 * t123;
t169 = pkin(4) * t123 - t198;
t168 = t124 * t171;
t165 = t123 * t129 + t198;
t164 = -t172 * MDP(26) - t98 * MDP(27);
t163 = t123 * t120 + t192;
t162 = -t123 * t137 - t192;
t161 = -t109 * MDP(26) - t99 * MDP(28);
t160 = -t150 * MDP(23) + t147 * MDP(24);
t159 = t150 * MDP(26) - t179;
t158 = -0.2e1 * t150 * MDP(28) - 0.2e1 * t178;
t157 = -t150 * t207 - t178;
t156 = (MDP(19) * t200 - t148 * MDP(20)) * pkin(3);
t155 = -t110 * MDP(20) + (t95 * t150 + t194) * MDP(29) + t208 * t109 + (-MDP(17) + t209) * t124 + ((t145 - t146) * MDP(22) - MDP(21) * t191 - MDP(16)) * t123;
t154 = (MDP(31) * qJ(6) + t206) * t150 + (-MDP(26) + t173) * t147;
t139 = t147 * MDP(29);
t1 = [MDP(1) + t121 * MDP(25) + (t95 ^ 2 + t96 ^ 2 + t99 ^ 2) * MDP(31) + (MDP(7) * t151 - 0.2e1 * t149 * MDP(8)) * t151 + ((-2 * MDP(4) + t199) * pkin(1)) + (t146 * MDP(21) + MDP(14) - 0.2e1 * t175) * t203 + (0.2e1 * t149 * MDP(12) + 0.2e1 * t151 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + 0.2e1 * (t134 * MDP(19) - t96 * MDP(28) + t95 * MDP(30) + t164) * t124 + 0.2e1 * (-t134 * MDP(20) + (MDP(15) + t160) * t124 + (t147 * t95 - t150 * t96) * MDP(29) + (-t147 * MDP(28) + t150 * MDP(30)) * t99 + (-t147 * MDP(26) - t150 * MDP(27)) * t109) * t123; MDP(4) - t199 + (t99 * t123 + t124 * t194) * MDP(31) + (t124 * t190 + t206 * t211) * t150 - t207 * t211 * t147; MDP(6) + (t121 * t184 + t203) * MDP(31); t99 * t181 + t155 + (MDP(26) * t162 - MDP(28) * t163 + t183 * t96 - t188) * t147 + (-t152 * MDP(13) - MDP(10)) * t149 + (MDP(27) * t162 + MDP(30) * t163 + t183 * t95 + t161) * t150 + (t152 * MDP(12) + MDP(9)) * t151; t151 * MDP(12) - t149 * MDP(13) + t136 * t168 + (t157 + t181) * t123 + t170; MDP(11) + t186 * t201 + t136 ^ 2 * t171 + (t158 + t181) * t120 + t174 - 0.2e1 * t159 * t137 + 0.2e1 * t156; t99 * t180 + (t169 * MDP(27) + t165 * MDP(30) + pkin(9) * t190 + t161) * t150 + (MDP(26) * t169 - MDP(28) * t165 + pkin(9) * t189 - t188) * t147 + t155; pkin(9) * t168 + (t157 + t180) * t123 + t170; (t185 + t186) * MDP(29) + (pkin(9) * t186 + t120 * t129) * MDP(31) + t156 + (MDP(26) * t196 + MDP(28) * t187) * t150 + (-MDP(27) * t196 + MDP(30) * t187) * t147 + t174; t185 * t201 + pkin(9) ^ 2 * t171 + (t158 + t180) * t129 + 0.2e1 * t159 * pkin(4) + t174; t124 * MDP(25) + (-t172 + 0.2e1 * t197) * MDP(28) + (t98 + 0.2e1 * t193) * MDP(30) + (-t96 * pkin(5) + t95 * qJ(6)) * MDP(31) + (MDP(29) * t167 + t160) * t123 + t164; t154 * t124; t136 * t154 + t176; pkin(9) * t154 + t176; MDP(25) + 0.2e1 * pkin(5) * MDP(28) + 0.2e1 * qJ(6) * MDP(30) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(31); -t150 * t123 * MDP(29) - t124 * MDP(28) + t189; t124 * t182; t136 * t182 + t139; pkin(9) * t182 + t139; t173; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
