% Calculate joint inertia matrix for
% S6RPRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR6_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:14:36
% EndTime: 2019-03-09 07:14:38
% DurationCPUTime: 0.72s
% Computational Cost: add. (1408->177), mult. (2711->245), div. (0->0), fcn. (3182->10), ass. (0->102)
t179 = sin(qJ(4));
t183 = cos(qJ(4));
t224 = pkin(8) + pkin(9);
t163 = t224 * t179;
t164 = t224 * t183;
t178 = sin(qJ(5));
t182 = cos(qJ(5));
t140 = -t182 * t163 - t164 * t178;
t141 = -t163 * t178 + t164 * t182;
t158 = t178 * t179 - t182 * t183;
t159 = t178 * t183 + t179 * t182;
t124 = -pkin(10) * t159 + t140;
t125 = -pkin(10) * t158 + t141;
t177 = sin(qJ(6));
t181 = cos(qJ(6));
t134 = t181 * t158 + t159 * t177;
t135 = -t158 * t177 + t159 * t181;
t194 = t135 * MDP(31) - t134 * MDP(32) + (t124 * t181 - t125 * t177) * MDP(34) - (t124 * t177 + t125 * t181) * MDP(35);
t186 = t159 * MDP(24) - t158 * MDP(25) + t140 * MDP(27) - t141 * MDP(28) + t194;
t190 = -MDP(20) * t179 - MDP(21) * t183;
t231 = t179 * MDP(17) + t183 * MDP(18) + pkin(8) * t190 + t186;
t149 = t158 * MDP(27);
t130 = t134 * MDP(34);
t207 = -t135 * MDP(35) - t130;
t193 = -t159 * MDP(28) - t149 + t207;
t176 = cos(pkin(11));
t167 = -pkin(2) * t176 - pkin(1);
t229 = 0.2e1 * t167;
t228 = -2 * MDP(23);
t227 = 0.2e1 * MDP(28);
t226 = -2 * MDP(30);
t225 = 0.2e1 * MDP(35);
t223 = cos(qJ(3));
t222 = pkin(1) * MDP(7);
t175 = sin(pkin(11));
t180 = sin(qJ(3));
t156 = t175 * t180 - t176 * t223;
t221 = pkin(4) * t156;
t220 = pkin(4) * t178;
t219 = pkin(4) * t182;
t218 = pkin(5) * t156;
t217 = pkin(5) * t181;
t216 = pkin(7) + qJ(2);
t157 = t175 * t223 + t180 * t176;
t129 = pkin(3) * t156 - pkin(8) * t157 + t167;
t161 = t216 * t175;
t162 = t216 * t176;
t137 = -t180 * t161 + t162 * t223;
t118 = t183 * t129 - t137 * t179;
t211 = t157 * t183;
t107 = -pkin(9) * t211 + t118 + t221;
t213 = t137 * t183;
t110 = t213 + (-pkin(9) * t157 + t129) * t179;
t214 = t110 * t182;
t105 = t107 * t178 + t214;
t126 = t159 * t157;
t103 = -pkin(10) * t126 + t105;
t215 = t103 * t181;
t212 = t157 * t179;
t210 = t175 * MDP(5);
t209 = t176 * MDP(4);
t208 = t179 * t183;
t205 = MDP(22) * t159;
t204 = MDP(29) * t135;
t169 = pkin(5) + t219;
t165 = t181 * t169;
t144 = -t177 * t220 + t165;
t203 = MDP(34) * t144;
t145 = t169 * t177 + t181 * t220;
t202 = MDP(35) * t145;
t201 = MDP(35) * t177;
t127 = t158 * t157;
t116 = t181 * t126 - t127 * t177;
t111 = t116 * MDP(32);
t117 = -t126 * t177 - t127 * t181;
t112 = t117 * MDP(31);
t122 = t126 * MDP(25);
t123 = t127 * MDP(24);
t200 = t157 * MDP(14);
t199 = t182 * MDP(27);
t198 = MDP(26) + MDP(33);
t197 = t156 * MDP(33) - t111 + t112;
t170 = -pkin(4) * t183 - pkin(3);
t196 = MDP(16) * t208;
t195 = MDP(19) + t198;
t104 = t182 * t107 - t110 * t178;
t102 = pkin(10) * t127 + t104 + t218;
t99 = t181 * t102 - t103 * t177;
t136 = t161 * t223 + t180 * t162;
t100 = t102 * t177 + t215;
t192 = MDP(17) * t183 - MDP(18) * t179;
t191 = MDP(20) * t183 - MDP(21) * t179;
t120 = pkin(4) * t212 + t136;
t189 = t156 * MDP(26) - t122 - t123 + t197;
t188 = -MDP(13) - t191;
t187 = (MDP(34) * t181 - t201) * pkin(5);
t174 = t183 ^ 2;
t173 = t179 ^ 2;
t142 = pkin(5) * t158 + t170;
t119 = t129 * t179 + t213;
t115 = t126 * pkin(5) + t120;
t1 = [t200 * t229 + MDP(1) - (-MDP(22) * t127 + t126 * t228) * t127 + (MDP(29) * t117 + t116 * t226) * t117 + (0.2e1 * t209 - 0.2e1 * t210 + t222) * pkin(1) + (t174 * MDP(15) + MDP(8) - 0.2e1 * t196) * t157 ^ 2 + t195 * t156 ^ 2 + (MDP(13) * t229 - 0.2e1 * t123 - 0.2e1 * t122 + 0.2e1 * t112 - 0.2e1 * t111 + 0.2e1 * (-MDP(9) + t192) * t157) * t156 + 0.2e1 * (t115 * t116 + t156 * t99) * MDP(34) + 0.2e1 * (t104 * t156 + t120 * t126) * MDP(27) + (-t100 * t156 + t115 * t117) * t225 + (-t105 * t156 - t120 * t127) * t227 + 0.2e1 * (t118 * t156 + t136 * t212) * MDP(20) + 0.2e1 * (-t119 * t156 + t136 * t211) * MDP(21) + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t175 ^ 2 + t176 ^ 2) * qJ(2); t200 - t209 + t210 - t222 + (-t188 + t193) * t156; MDP(7); -t137 * MDP(14) - t127 * t205 + (-t126 * t159 + t127 * t158) * MDP(23) + (t120 * t158 + t126 * t170) * MDP(27) + (t120 * t159 - t127 * t170) * MDP(28) + t117 * t204 + (-t116 * t135 - t117 * t134) * MDP(30) + (t115 * t134 + t116 * t142) * MDP(34) + (t115 * t135 + t117 * t142) * MDP(35) + t188 * t136 + (MDP(10) + MDP(15) * t208 + (-t173 + t174) * MDP(16) + t190 * pkin(3)) * t157 + (-MDP(11) + t231) * t156; 0; 0.2e1 * t196 + 0.2e1 * t170 * t149 + 0.2e1 * t142 * t130 + MDP(15) * t173 + MDP(12) + 0.2e1 * t191 * pkin(3) + (t158 * t228 + t170 * t227 + t205) * t159 + (t134 * t226 + t142 * t225 + t204) * t135; t156 * MDP(19) + t118 * MDP(20) - t119 * MDP(21) + (t156 * t219 + t104) * MDP(27) + (-t214 + (-t107 - t221) * t178) * MDP(28) + (t144 * t156 + t99) * MDP(34) + (-t145 * t156 - t100) * MDP(35) + t192 * t157 + t189; t191 + t193; t231; 0.2e1 * (-MDP(28) * t178 + t199) * pkin(4) + 0.2e1 * t203 - 0.2e1 * t202 + t195; t104 * MDP(27) - t105 * MDP(28) + (t156 * t217 + t99) * MDP(34) + (-t215 + (-t102 - t218) * t177) * MDP(35) + t189; t193; t186; (t165 + t217) * MDP(34) + (-pkin(5) - t169) * t201 + (t199 + (-MDP(34) * t177 - MDP(35) * t181 - MDP(28)) * t178) * pkin(4) + t198; 0.2e1 * t187 + t198; t99 * MDP(34) - t100 * MDP(35) + t197; t207; t194; MDP(33) - t202 + t203; MDP(33) + t187; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
