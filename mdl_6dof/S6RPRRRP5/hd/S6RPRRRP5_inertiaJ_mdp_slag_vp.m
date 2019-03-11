% Calculate joint inertia matrix for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRP5_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:12:27
% EndTime: 2019-03-09 06:12:30
% DurationCPUTime: 0.73s
% Computational Cost: add. (1430->164), mult. (2605->219), div. (0->0), fcn. (3015->8), ass. (0->88)
t153 = sin(qJ(4));
t138 = t153 * pkin(3) + pkin(9);
t152 = sin(qJ(5));
t148 = t152 ^ 2;
t154 = cos(qJ(5));
t149 = t154 ^ 2;
t192 = t148 + t149;
t195 = t192 * t138;
t216 = t152 * MDP(24) + t154 * MDP(25);
t150 = sin(pkin(10));
t151 = cos(pkin(10));
t208 = sin(qJ(3));
t210 = cos(qJ(3));
t127 = t210 * t150 + t208 * t151;
t176 = t208 * t150 - t210 * t151;
t209 = cos(qJ(4));
t122 = t209 * t127 - t153 * t176;
t215 = 0.2e1 * t122;
t214 = t152 * MDP(28);
t137 = -t151 * pkin(2) - pkin(1);
t123 = t176 * pkin(3) + t137;
t213 = 0.2e1 * t123;
t212 = 0.2e1 * t137;
t211 = 2 * MDP(30);
t207 = pkin(1) * MDP(7);
t121 = t153 * t127 + t209 * t176;
t206 = pkin(9) * t121;
t205 = t121 * pkin(5);
t184 = t209 * pkin(3);
t139 = -t184 - pkin(4);
t204 = pkin(4) - t139;
t203 = pkin(7) + qJ(2);
t202 = t121 * qJ(6);
t201 = t121 * t138;
t200 = t150 * MDP(5);
t199 = t151 * MDP(4);
t198 = t152 * t154;
t130 = t203 * t150;
t131 = t203 * t151;
t175 = -t210 * t130 - t208 * t131;
t113 = -t127 * pkin(8) + t175;
t160 = t208 * t130 - t210 * t131;
t114 = -t176 * pkin(8) - t160;
t110 = t153 * t113 + t209 * t114;
t111 = t121 * pkin(4) - t122 * pkin(9) + t123;
t179 = t152 * t110 - t154 * t111;
t99 = t179 - t205;
t197 = t99 * MDP(32);
t101 = t154 * t110 + t152 * t111;
t173 = -t154 * pkin(5) - t152 * qJ(6);
t129 = -pkin(4) + t173;
t125 = -t184 + t129;
t196 = -t125 - t129;
t194 = t192 * pkin(9);
t191 = MDP(32) * t129;
t190 = MDP(32) * t138;
t189 = MDP(32) * t152;
t109 = -t209 * t113 + t153 * t114;
t172 = pkin(5) * t152 - t154 * qJ(6);
t102 = t172 * t122 + t109;
t188 = t102 * MDP(31);
t187 = t121 * MDP(26);
t186 = t125 * MDP(32);
t185 = -MDP(28) + MDP(31);
t183 = -t172 * MDP(30) + t216;
t182 = MDP(23) * t198;
t181 = t148 * MDP(22) + MDP(19) + 0.2e1 * t182;
t180 = -MDP(32) * pkin(5) - MDP(29);
t178 = t192 * MDP(32);
t177 = -pkin(4) * t122 - t206;
t174 = MDP(27) - t180;
t171 = -t122 * t129 + t206;
t98 = t101 + t202;
t170 = t98 * t152 - t99 * t154;
t169 = MDP(32) * qJ(6) + t185;
t168 = -t122 * t125 + t201;
t167 = t122 * t139 - t201;
t166 = t154 * MDP(24) - t152 * MDP(25);
t165 = -t179 * MDP(27) - t101 * MDP(28);
t164 = t154 * MDP(27) - t214;
t163 = -t109 * MDP(27) - t102 * MDP(29);
t162 = -0.2e1 * t154 * MDP(29) - 0.2e1 * t152 * MDP(31);
t161 = t176 * MDP(13);
t159 = (t209 * MDP(20) - t153 * MDP(21)) * pkin(3);
t158 = -t110 * MDP(21) + (t99 * t152 + t98 * t154) * MDP(30) + (t214 - MDP(20)) * t109 + ((-t148 + t149) * MDP(23) + MDP(22) * t198 + MDP(17)) * t122 + (-MDP(18) + t216) * t121;
t157 = -t174 * t152 + t169 * t154;
t141 = t152 * MDP(30);
t1 = [(t102 ^ 2 + t98 ^ 2 + t99 ^ 2) * MDP(32) + t161 * t212 + MDP(1) + (0.2e1 * t199 - 0.2e1 * t200 + t207) * pkin(1) + (MDP(14) * t212 + MDP(8) * t127 - 0.2e1 * t176 * MDP(9)) * t127 + (MDP(20) * t213 + t187 + (-MDP(16) + t166) * t215) * t121 + 0.2e1 * (-t99 * MDP(29) + t98 * MDP(31) + t165) * t121 + (-t170 * MDP(30) + (t152 * MDP(27) + t154 * MDP(28)) * t109 + (t152 * MDP(29) - t154 * MDP(31)) * t102) * t215 + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t150 ^ 2 + t151 ^ 2) * qJ(2) + (MDP(21) * t213 + (t149 * MDP(22) + MDP(15) - 0.2e1 * t182) * t122) * t122; -t199 + t200 - t207 + t161 + t127 * MDP(14) + t170 * MDP(32) + (-t192 * MDP(30) + MDP(21)) * t122 + (MDP(20) + (MDP(27) + MDP(29)) * t154 + t185 * t152) * t121; MDP(7) + t178; t158 + t102 * t186 + (t167 * MDP(28) + t168 * MDP(31) + t98 * t190 + t163) * t154 + (t167 * MDP(27) - t168 * MDP(29) + t99 * t190 - t188) * t152 - t176 * MDP(11) + t127 * MDP(10) + t175 * MDP(13) + t160 * MDP(14); 0; MDP(12) + t195 * t211 + t138 ^ 2 * t178 + (t162 + t186) * t125 + t181 - 0.2e1 * t164 * t139 + 0.2e1 * t159; t102 * t191 + (pkin(9) * t98 * MDP(32) + t177 * MDP(28) + t171 * MDP(31) + t163) * t154 + (t177 * MDP(27) - t171 * MDP(29) + pkin(9) * t197 - t188) * t152 + t158; 0; (t194 + t195) * MDP(30) + (pkin(9) * t195 + t125 * t129) * MDP(32) + t159 + (t204 * MDP(27) + t196 * MDP(29)) * t154 + (-t204 * MDP(28) + t196 * MDP(31)) * t152 + t181; t194 * t211 + pkin(9) ^ 2 * t178 + (t162 + t191) * t129 + 0.2e1 * t164 * pkin(4) + t181; t187 + (-t179 + 0.2e1 * t205) * MDP(29) + (t101 + 0.2e1 * t202) * MDP(31) + (-t99 * pkin(5) + t98 * qJ(6)) * MDP(32) + (t173 * MDP(30) + t166) * t122 + t165; t169 * t152 + t174 * t154; t157 * t138 + t183; t157 * pkin(9) + t183; MDP(26) + 0.2e1 * pkin(5) * MDP(29) + 0.2e1 * qJ(6) * MDP(31) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(32); t154 * t122 * MDP(30) - t121 * MDP(29) + t197; -t154 * MDP(32); t138 * t189 + t141; pkin(9) * t189 + t141; t180; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
