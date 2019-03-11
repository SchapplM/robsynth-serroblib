% Calculate joint inertia matrix for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRRP11_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:41:16
% EndTime: 2019-03-09 06:41:20
% DurationCPUTime: 1.25s
% Computational Cost: add. (2574->249), mult. (6792->379), div. (0->0), fcn. (7731->12), ass. (0->106)
t171 = sin(pkin(6));
t170 = sin(pkin(7));
t214 = cos(pkin(6));
t189 = t214 * t170;
t172 = cos(pkin(12));
t173 = cos(pkin(7));
t206 = t172 * t173;
t228 = t171 * t206 + t189;
t174 = sin(qJ(5));
t177 = cos(qJ(5));
t183 = t174 * MDP(27) + t177 * MDP(28);
t227 = t174 * MDP(24) + t177 * MDP(25) - pkin(11) * t183 - MDP(18);
t169 = sin(pkin(12));
t191 = pkin(1) * t214;
t207 = t171 * t172;
t149 = qJ(2) * t207 + t169 * t191;
t133 = t228 * pkin(9) + t149;
t176 = sin(qJ(3));
t179 = cos(qJ(3));
t161 = t172 * t191;
t210 = t169 * t171;
t137 = t214 * pkin(2) + t161 + (-pkin(9) * t173 - qJ(2)) * t210;
t142 = (-pkin(9) * t169 * t170 - pkin(2) * t172 - pkin(1)) * t171;
t186 = t137 * t173 + t142 * t170;
t124 = -t176 * t133 + t179 * t186;
t226 = 0.2e1 * MDP(27);
t225 = 0.2e1 * MDP(28);
t224 = 2 * MDP(29);
t165 = t171 ^ 2;
t223 = pkin(1) * t165;
t222 = pkin(10) * t174;
t221 = pkin(10) * t177;
t178 = cos(qJ(4));
t220 = pkin(10) * t178;
t219 = -qJ(6) - pkin(11);
t218 = MDP(30) * pkin(5);
t217 = pkin(3) * MDP(20);
t216 = pkin(3) * MDP(21);
t215 = pkin(10) * MDP(21);
t175 = sin(qJ(4));
t213 = qJ(6) * t175;
t128 = -t170 * t137 + t173 * t142;
t135 = t176 * t210 - t228 * t179;
t136 = t176 * t189 + (t169 * t179 + t176 * t206) * t171;
t120 = t135 * pkin(3) - t136 * pkin(10) + t128;
t125 = t179 * t133 + t176 * t186;
t146 = t170 * t207 - t214 * t173;
t123 = -t146 * pkin(10) + t125;
t116 = t178 * t120 - t175 * t123;
t114 = -t135 * pkin(4) - t116;
t212 = t114 * t174;
t211 = t114 * t177;
t209 = t170 * t176;
t208 = t170 * t179;
t205 = t174 * t177;
t130 = t136 * t178 - t146 * t175;
t126 = t130 * t174 - t135 * t177;
t204 = t126 * MDP(25);
t127 = t130 * t177 + t135 * t174;
t203 = t127 * MDP(22);
t202 = t127 * MDP(24);
t129 = t136 * t175 + t146 * t178;
t201 = t129 * MDP(26);
t200 = t130 * MDP(16);
t199 = t130 * MDP(17);
t198 = t135 * MDP(19);
t151 = t175 * t173 + t178 * t209;
t139 = t177 * t151 - t174 * t208;
t197 = t139 * MDP(28);
t164 = -t177 * pkin(5) - pkin(4);
t196 = t164 * MDP(30);
t195 = t174 * MDP(25);
t194 = t178 * MDP(26);
t193 = t177 * t220;
t190 = MDP(23) * t205;
t188 = -pkin(10) * MDP(20) + MDP(17);
t187 = -MDP(29) * pkin(5) + MDP(24);
t117 = t175 * t120 + t178 * t123;
t115 = t135 * pkin(11) + t117;
t122 = t146 * pkin(3) - t124;
t119 = t129 * pkin(4) - t130 * pkin(11) + t122;
t111 = -t174 * t115 + t177 * t119;
t112 = t177 * t115 + t174 * t119;
t155 = -t178 * pkin(4) - t175 * pkin(11) - pkin(3);
t152 = t177 * t155;
t143 = -t174 * t220 + t152;
t144 = t174 * t155 + t193;
t185 = t143 * MDP(27) - t144 * MDP(28);
t184 = t177 * MDP(27) - t174 * MDP(28);
t182 = t177 * MDP(24) - MDP(16) - t195;
t180 = t111 * MDP(27) - t112 * MDP(28) + t201 + t202 - t204;
t168 = t177 ^ 2;
t167 = t175 ^ 2;
t166 = t174 ^ 2;
t157 = t219 * t177;
t156 = t219 * t174;
t154 = (pkin(5) * t174 + pkin(10)) * t175;
t150 = -t178 * t173 + t175 * t209;
t148 = -qJ(2) * t210 + t161;
t140 = t193 + (t155 - t213) * t174;
t138 = -t174 * t151 - t177 * t208;
t134 = -t177 * t213 + t152 + (-pkin(5) - t222) * t178;
t113 = t126 * pkin(5) + t114;
t110 = -t126 * qJ(6) + t112;
t109 = t129 * pkin(5) - t127 * qJ(6) + t111;
t1 = [t146 ^ 2 * MDP(12) + MDP(1) + (t109 ^ 2 + t110 ^ 2 + t113 ^ 2) * MDP(30) + t130 ^ 2 * MDP(15) + (t165 * pkin(1) ^ 2 + t148 ^ 2 + t149 ^ 2) * MDP(7) + (-0.2e1 * t146 * MDP(10) + MDP(8) * t136) * t136 + (-0.2e1 * t126 * MDP(23) + t203) * t127 + (0.2e1 * t146 * MDP(11) - 0.2e1 * t136 * MDP(9) + t198 + 0.2e1 * t199) * t135 + (-0.2e1 * t135 * MDP(18) - 0.2e1 * t200 + t201 + 0.2e1 * t202 - 0.2e1 * t204) * t129 + 0.2e1 * (t125 * t146 + t128 * t136) * MDP(14) + 0.2e1 * (-t124 * t146 + t128 * t135) * MDP(13) + 0.2e1 * (t148 * t214 + t172 * t223) * MDP(4) + 0.2e1 * (-t149 * t214 - t169 * t223) * MDP(5) + (-t112 * t129 + t114 * t127) * t225 + 0.2e1 * (t116 * t135 + t122 * t129) * MDP(20) + 0.2e1 * (-t117 * t135 + t122 * t130) * MDP(21) + (-t109 * t127 - t110 * t126) * t224 + (t111 * t129 + t114 * t126) * t226 + 0.2e1 * (-t148 * t169 + t149 * t172) * MDP(6) * t171; (t173 * t135 - t146 * t208) * MDP(13) + (t173 * t136 + t146 * t209) * MDP(14) + (-t129 * t208 - t150 * t135) * MDP(20) + (-t130 * t208 - t151 * t135) * MDP(21) + (t150 * t126 + t138 * t129) * MDP(27) + (t150 * t127 - t139 * t129) * MDP(28) + (-t139 * t126 - t138 * t127) * MDP(29) + (t109 * t138 + t110 * t139 + t113 * t150) * MDP(30) + (-t172 * MDP(4) + t169 * MDP(5) - pkin(1) * MDP(7)) * t171; MDP(7) + (t138 ^ 2 + t139 ^ 2 + t150 ^ 2) * MDP(30); t136 * MDP(10) - t135 * MDP(11) - t146 * MDP(12) + t124 * MDP(13) - t125 * MDP(14) - t130 * t216 + (-t140 * t126 - t134 * t127) * MDP(29) + (t109 * t134 + t110 * t140 + t113 * t154) * MDP(30) + (t185 - t217) * t129 + (t200 - t122 * MDP(20) + (MDP(18) - t215) * t135 - t180) * t178 + (t130 * MDP(15) + t122 * MDP(21) + t177 * t203 + (-t126 * t177 - t127 * t174) * MDP(23) + (pkin(10) * t126 + t212) * MDP(27) + (pkin(10) * t127 + t211) * MDP(28) + (-t109 * t177 - t110 * t174) * MDP(29) + t188 * t135 + t182 * t129) * t175; (t138 * t134 + t139 * t140 + t150 * t154) * MDP(30) + (-t138 * MDP(27) + t197) * t178 + ((-t138 * t177 - t139 * t174) * MDP(29) + t183 * t150) * t175 + (-t176 * MDP(14) + (MDP(20) * t178 - MDP(21) * t175 + MDP(13)) * t179) * t170; MDP(12) + (t134 ^ 2 + t140 ^ 2 + t154 ^ 2) * MDP(30) + (t194 + 0.2e1 * t217) * t178 + (t168 * MDP(22) + MDP(15) - 0.2e1 * t190) * t167 + (-t143 * t178 + t167 * t222) * t226 + (t144 * t178 + t167 * t221) * t225 + (-0.2e1 * t216 + (-t134 * t177 - t140 * t174) * t224 - 0.2e1 * t182 * t178) * t175; t199 + t198 + t116 * MDP(20) - t117 * MDP(21) + t174 * t203 + (-t174 * t126 + t127 * t177) * MDP(23) + (-pkin(4) * t126 - t211) * MDP(27) + (-pkin(4) * t127 + t212) * MDP(28) + (-t109 * t174 + t110 * t177 + t157 * t126 - t156 * t127) * MDP(29) + (t109 * t156 - t110 * t157 + t113 * t164) * MDP(30) + t227 * t129; -t151 * MDP(21) + (-t138 * t174 + t139 * t177) * MDP(29) + (t138 * t156 - t139 * t157) * MDP(30) + (-MDP(20) - t184 + t196) * t150; (-t134 * t174 + t140 * t177) * MDP(29) + (t134 * t156 - t140 * t157 + t154 * t164) * MDP(30) + (-t215 - t227) * t178 + (MDP(22) * t205 + (-t166 + t168) * MDP(23) + (-pkin(4) * t174 - t221) * MDP(27) + (-pkin(4) * t177 + t222) * MDP(28) + (-t156 * t177 + t157 * t174) * MDP(29) + t188) * t175; MDP(19) + t166 * MDP(22) + 0.2e1 * t190 + (-t156 * t174 - t157 * t177) * t224 + (t156 ^ 2 + t157 ^ 2 + t164 ^ 2) * MDP(30) + 0.2e1 * t184 * pkin(4); (-t127 * MDP(29) + t109 * MDP(30)) * pkin(5) + t180; -t197 + (MDP(27) + t218) * t138; t134 * t218 - t194 + (t187 * t177 - t195) * t175 + t185; t156 * t218 + (-MDP(28) * pkin(11) + MDP(25)) * t177 + (-MDP(27) * pkin(11) + t187) * t174; MDP(30) * pkin(5) ^ 2 + MDP(26); t113 * MDP(30); t150 * MDP(30); t154 * MDP(30); t196; 0; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
