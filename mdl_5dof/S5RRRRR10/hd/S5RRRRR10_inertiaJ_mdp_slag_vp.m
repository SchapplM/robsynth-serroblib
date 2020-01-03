% Calculate joint inertia matrix for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR10_inertiaJ_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:36:00
% EndTime: 2019-12-31 22:36:03
% DurationCPUTime: 0.84s
% Computational Cost: add. (1054->193), mult. (2401->280), div. (0->0), fcn. (2654->10), ass. (0->98)
t166 = sin(qJ(5));
t170 = cos(qJ(5));
t203 = t166 * MDP(27) + t170 * MDP(28);
t165 = cos(pkin(5));
t168 = sin(qJ(3));
t171 = cos(qJ(3));
t164 = sin(pkin(5));
t169 = sin(qJ(2));
t207 = t164 * t169;
t140 = -t165 * t171 + t168 * t207;
t141 = t165 * t168 + t171 * t207;
t167 = sin(qJ(4));
t216 = cos(qJ(4));
t122 = t216 * t140 + t141 * t167;
t123 = -t167 * t140 + t216 * t141;
t225 = t123 * MDP(20) - t122 * MDP(21);
t179 = MDP(30) * t170 - MDP(31) * t166;
t224 = -(MDP(16) * t168 + MDP(17) * t171) * pkin(8) + t168 * MDP(13) + t171 * MDP(14);
t217 = pkin(8) + pkin(9);
t150 = t217 * t171;
t188 = t216 * t168;
t130 = t150 * t167 + t217 * t188;
t205 = t167 * t168;
t131 = t216 * t150 - t217 * t205;
t147 = -t216 * t171 + t205;
t148 = t167 * t171 + t188;
t223 = t148 * MDP(20) - t147 * MDP(21) - t130 * MDP(23) - t131 * MDP(24);
t222 = 0.2e1 * MDP(16);
t221 = 0.2e1 * MDP(23);
t220 = 0.2e1 * MDP(24);
t219 = 0.2e1 * MDP(30);
t218 = 0.2e1 * MDP(31);
t215 = pkin(1) * t169;
t172 = cos(qJ(2));
t214 = pkin(1) * t172;
t206 = t164 * t172;
t191 = pkin(7) * t206;
t135 = t191 + (pkin(8) + t215) * t165;
t136 = (-pkin(2) * t172 - pkin(8) * t169 - pkin(1)) * t164;
t116 = -t168 * t135 + t171 * t136;
t192 = pkin(3) * t206;
t109 = -t141 * pkin(9) + t116 - t192;
t117 = t135 * t171 + t136 * t168;
t112 = -pkin(9) * t140 + t117;
t104 = t216 * t109 - t167 * t112;
t102 = pkin(4) * t206 - t104;
t212 = t102 * t170;
t211 = t122 * t166;
t210 = t122 * t170;
t209 = t130 * t170;
t208 = t148 * t166;
t204 = t169 * MDP(6);
t202 = MDP(29) * t147;
t114 = t166 * t123 + t170 * t206;
t199 = t114 * MDP(28);
t115 = t170 * t123 - t166 * t206;
t198 = t115 * MDP(25);
t197 = t115 * MDP(27);
t196 = t122 * MDP(29);
t195 = t123 * MDP(19);
t194 = t141 * MDP(11);
t193 = t170 * MDP(25);
t190 = t216 * pkin(3);
t156 = -pkin(3) * t171 - pkin(2);
t189 = t216 * t112;
t187 = t166 * t170 * MDP(26);
t161 = t166 ^ 2;
t186 = t161 * MDP(25) + MDP(22) + 0.2e1 * t187;
t185 = -pkin(4) * t148 - pkin(10) * t147;
t154 = pkin(3) * t167 + pkin(10);
t155 = -t190 - pkin(4);
t184 = -t147 * t154 + t148 * t155;
t183 = t141 * MDP(13) - t140 * MDP(14);
t180 = MDP(27) * t170 - MDP(28) * t166;
t178 = t166 * MDP(30) + t170 * MDP(31);
t152 = pkin(7) * t207;
t134 = t152 + (-pkin(2) - t214) * t165;
t105 = t167 * t109 + t189;
t177 = -MDP(19) + t180;
t176 = (t216 * MDP(23) - t167 * MDP(24)) * pkin(3);
t162 = t170 ^ 2;
t175 = (-t161 + t162) * t148 * MDP(26) + t193 * t208 + t223 + t203 * t147;
t124 = t140 * pkin(3) + t134;
t174 = -MDP(22) * t206 + (-t114 * t166 + t115 * t170) * MDP(26) + t166 * t198 + t225 + t203 * t122;
t103 = -pkin(10) * t206 + t105;
t106 = t122 * pkin(4) - t123 * pkin(10) + t124;
t100 = t103 * t170 + t106 * t166;
t99 = -t103 * t166 + t106 * t170;
t173 = t99 * MDP(30) - t100 * MDP(31) + t196 + t197 - t199;
t160 = t164 ^ 2;
t143 = t165 * t215 + t191;
t142 = t165 * t214 - t152;
t127 = t130 * t166;
t126 = pkin(4) * t147 - pkin(10) * t148 + t156;
t111 = t126 * t166 + t131 * t170;
t110 = t126 * t170 - t131 * t166;
t101 = t102 * t166;
t1 = [t123 ^ 2 * MDP(18) + t165 ^ 2 * MDP(8) + MDP(1) + (-0.2e1 * t140 * MDP(12) + t194) * t141 + (-0.2e1 * t114 * MDP(26) + t198) * t115 + ((MDP(4) * t169 + 0.2e1 * MDP(5) * t172) * t169 + (MDP(15) + MDP(22)) * t172 ^ 2) * t160 + (-0.2e1 * t195 + t196 + 0.2e1 * t197 - 0.2e1 * t199) * t122 + 0.2e1 * (t165 * t204 + (MDP(7) * t165 - t183 - t225) * t172) * t164 + (t102 * t114 + t122 * t99) * t219 + (-t100 * t122 + t102 * t115) * t218 + 0.2e1 * (-t143 * t165 - t160 * t215) * MDP(10) + (-t104 * t206 + t124 * t122) * t221 + (t105 * t206 + t124 * t123) * t220 + (-t116 * t206 + t134 * t140) * t222 + 0.2e1 * (t117 * t206 + t134 * t141) * MDP(17) + 0.2e1 * (t142 * t165 + t160 * t214) * MDP(9); t165 * MDP(8) + t142 * MDP(9) - t143 * MDP(10) + (-t140 * t168 + t141 * t171) * MDP(12) + (-pkin(2) * t140 - t134 * t171) * MDP(16) + (-pkin(2) * t141 + t134 * t168) * MDP(17) + (t110 * t122 + t114 * t130) * MDP(30) + (-t111 * t122 + t115 * t130) * MDP(31) + t168 * t194 + (t122 * MDP(23) + t123 * MDP(24)) * t156 + (t124 * MDP(23) + t173 - t195) * t147 + (t124 * MDP(24) + (-t114 * t170 - t115 * t166) * MDP(26) + t123 * MDP(18) + t115 * t193 + t178 * t102 + t177 * t122) * t148 + (t204 + (MDP(7) - t223 - t224) * t172) * t164; pkin(2) * t171 * t222 + MDP(8) + t130 * t208 * t219 + (MDP(11) * t168 + 0.2e1 * t171 * MDP(12) - 0.2e1 * pkin(2) * MDP(17)) * t168 + (t110 * t219 - t111 * t218 + t156 * t221 + t202) * t147 + (0.2e1 * t177 * t147 + t156 * t220 + t209 * t218 + (t162 * MDP(25) + MDP(18) - 0.2e1 * t187) * t148) * t148; -MDP(15) * t206 + t116 * MDP(16) - t117 * MDP(17) + (-t190 * t206 + t104) * MDP(23) + (-t189 + (-t109 + t192) * t167) * MDP(24) + (t114 * t155 - t154 * t211 - t212) * MDP(30) + (t115 * t155 - t154 * t210 + t101) * MDP(31) + t174 + t183; (t166 * t184 - t209) * MDP(30) + (t170 * t184 + t127) * MDP(31) + t175 + t224; -0.2e1 * t155 * t179 + MDP(15) + 0.2e1 * t176 + t186; t104 * MDP(23) - t105 * MDP(24) + (-pkin(4) * t114 - pkin(10) * t211 - t212) * MDP(30) + (-pkin(4) * t115 - pkin(10) * t210 + t101) * MDP(31) + t174; (t166 * t185 - t209) * MDP(30) + (t170 * t185 + t127) * MDP(31) + t175; t176 + t186 + t179 * (pkin(4) - t155); 0.2e1 * pkin(4) * t179 + t186; t173; t110 * MDP(30) - t111 * MDP(31) + t148 * t180 + t202; -t154 * t178 + t203; -pkin(10) * t178 + t203; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
