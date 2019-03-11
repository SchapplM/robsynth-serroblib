% Calculate joint inertia matrix for
% S6RRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP1_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:58:16
% EndTime: 2019-03-10 00:58:18
% DurationCPUTime: 0.68s
% Computational Cost: add. (1451->199), mult. (2594->256), div. (0->0), fcn. (2920->8), ass. (0->87)
t162 = sin(qJ(5));
t166 = cos(qJ(5));
t195 = t162 * MDP(27) + t166 * MDP(28);
t186 = t166 * MDP(30);
t207 = t162 * MDP(31);
t172 = 0.2e1 * t186 - 0.2e1 * t207;
t203 = 2 * MDP(32);
t164 = sin(qJ(3));
t165 = sin(qJ(2));
t168 = cos(qJ(3));
t169 = cos(qJ(2));
t138 = t164 * t165 - t168 * t169;
t154 = -t169 * pkin(2) - pkin(1);
t126 = t138 * pkin(3) + t154;
t206 = 0.2e1 * t126;
t205 = 0.2e1 * t154;
t204 = 0.2e1 * t169;
t202 = pkin(7) + pkin(8);
t201 = pkin(2) * t164;
t200 = pkin(4) * t162;
t167 = cos(qJ(4));
t199 = t167 * pkin(3);
t198 = t162 * t166;
t158 = t166 * qJ(6);
t139 = t164 * t169 + t168 * t165;
t144 = t202 * t165;
t145 = t202 * t169;
t182 = -t168 * t144 - t164 * t145;
t110 = -t139 * pkin(9) + t182;
t176 = t164 * t144 - t168 * t145;
t111 = -t138 * pkin(9) - t176;
t163 = sin(qJ(4));
t107 = t163 * t110 + t167 * t111;
t197 = t166 * t107;
t152 = t168 * pkin(2) + pkin(3);
t196 = -t167 * t152 + t163 * t201;
t194 = MDP(28) * t162;
t118 = -t163 * t138 + t167 * t139;
t193 = MDP(32) * t118;
t106 = -t167 * t110 + t163 * t111;
t192 = t106 * MDP(30);
t117 = t167 * t138 + t163 * t139;
t191 = t117 * MDP(29);
t190 = t196 * MDP(23);
t131 = -t163 * t152 - t167 * t201;
t189 = t131 * MDP(24);
t188 = t162 * MDP(32);
t187 = t163 * MDP(24);
t185 = t168 * MDP(16);
t153 = -t166 * pkin(5) - pkin(4);
t184 = MDP(26) * t198;
t160 = t162 ^ 2;
t183 = t160 * MDP(25) + MDP(22) + 0.2e1 * t184;
t105 = t117 * pkin(4) - t118 * pkin(10) + t126;
t98 = t166 * t105 - t162 * t107;
t181 = MDP(15) + t183;
t180 = -pkin(4) * t118 - pkin(10) * t117;
t179 = t98 * MDP(30) - (t162 * t105 + t197) * MDP(31);
t128 = -pkin(4) + t196;
t129 = pkin(10) - t131;
t178 = -t117 * t129 + t118 * t128;
t150 = t163 * pkin(3) + pkin(10);
t151 = -pkin(4) - t199;
t177 = -t117 * t150 + t118 * t151;
t174 = t162 * MDP(30) + t166 * MDP(31);
t173 = (t167 * MDP(23) - t187) * pkin(3);
t161 = t166 ^ 2;
t97 = t197 + (-qJ(6) * t118 + t105) * t162;
t171 = t97 * t166 * MDP(32) - t107 * MDP(24) + (t207 - MDP(23)) * t106 + ((-t160 + t161) * MDP(26) + MDP(25) * t198 + MDP(20)) * t118 + (-MDP(21) + t195) * t117;
t170 = t139 * MDP(13) - t138 * MDP(14) + t182 * MDP(16) + t176 * MDP(17) + t171;
t159 = pkin(4) * t166;
t147 = t151 * t162;
t143 = t166 * pkin(10) + t158;
t142 = (-qJ(6) - pkin(10)) * t162;
t141 = t153 - t199;
t137 = t143 * t166;
t136 = t166 * t150 + t158;
t135 = (-qJ(6) - t150) * t162;
t127 = t136 * t166;
t125 = t128 * t162;
t124 = t153 + t196;
t123 = t166 * t129 + t158;
t122 = (-qJ(6) - t129) * t162;
t119 = t123 * t166;
t100 = t162 * t118 * pkin(5) + t106;
t95 = t117 * pkin(5) - t118 * t158 + t98;
t1 = [MDP(1) + pkin(1) * MDP(9) * t204 + t138 * MDP(16) * t205 + (t100 ^ 2 + t95 ^ 2 + t97 ^ 2) * MDP(33) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t165 + MDP(5) * t204) * t165 + (MDP(11) * t139 - 0.2e1 * t138 * MDP(12) + MDP(17) * t205) * t139 + (MDP(24) * t206 + (-t162 * t97 - t166 * t95) * t203 + 0.2e1 * t174 * t106 + (t161 * MDP(25) + MDP(18) - 0.2e1 * t184) * t118) * t118 + (MDP(23) * t206 + t191 + 0.2e1 * (MDP(27) * t166 - MDP(19) - t194) * t118 + 0.2e1 * t179) * t117; (t178 * MDP(30) + (-t118 * t123 - t95) * MDP(32)) * t162 + (MDP(31) * t178 - t122 * t193 - t192) * t166 + t170 + (-t169 * MDP(10) - t165 * MDP(9)) * pkin(7) + (t100 * t124 + t95 * t122 + t97 * t123) * MDP(33) + t165 * MDP(6) + t169 * MDP(7); MDP(8) + (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) * MDP(33) - t128 * t172 + 0.2e1 * (-t164 * MDP(17) + t185) * pkin(2) - 0.2e1 * t190 + 0.2e1 * t189 + (-t122 * t162 + t119) * t203 + t181; t170 + (MDP(31) * t177 - t135 * t193 - t192) * t166 + (t177 * MDP(30) + (-t118 * t136 - t95) * MDP(32)) * t162 + (t100 * t141 + t95 * t135 + t97 * t136) * MDP(33); (-t196 + t199) * MDP(23) + (t147 + t125) * MDP(31) + (t119 + t127 + (-t122 - t135) * t162) * MDP(32) + (t122 * t135 + t123 * t136 + t124 * t141) * MDP(33) + (-t128 - t151) * t186 + (-pkin(3) - t152) * t187 + (t185 + (-MDP(24) * t167 - MDP(17)) * t164) * pkin(2) + t181; (-t135 * t162 + t127) * t203 + (t135 ^ 2 + t136 ^ 2 + t141 ^ 2) * MDP(33) - t151 * t172 + 0.2e1 * t173 + t181; (t100 * t153 + t95 * t142 + t97 * t143) * MDP(33) + (MDP(31) * t180 - t142 * t193 - t192) * t166 + (t180 * MDP(30) + (-t118 * t143 - t95) * MDP(32)) * t162 + t171; -t190 + t189 + (-t128 * t166 + t159) * MDP(30) + (t125 - t200) * MDP(31) + (t119 + t137 + (-t122 - t142) * t162) * MDP(32) + (t122 * t142 + t123 * t143 + t124 * t153) * MDP(33) + t183; (-t151 * t166 + t159) * MDP(30) + (t147 - t200) * MDP(31) + (t127 + t137 + (-t135 - t142) * t162) * MDP(32) + (t135 * t142 + t136 * t143 + t141 * t153) * MDP(33) + t173 + t183; (-t142 * t162 + t137) * t203 + (t142 ^ 2 + t143 ^ 2 + t153 ^ 2) * MDP(33) + pkin(4) * t172 + t183; t95 * pkin(5) * MDP(33) + t191 + (-t194 + (-MDP(32) * pkin(5) + MDP(27)) * t166) * t118 + t179; -t174 * t129 + (t122 * MDP(33) - t188) * pkin(5) + t195; -t174 * t150 + (t135 * MDP(33) - t188) * pkin(5) + t195; -t174 * pkin(10) + (MDP(33) * t142 - t188) * pkin(5) + t195; MDP(33) * pkin(5) ^ 2 + MDP(29); t100 * MDP(33); t124 * MDP(33); t141 * MDP(33); t153 * MDP(33); 0; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
