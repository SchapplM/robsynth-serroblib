% Calculate joint inertia matrix for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPPRR3_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:58:40
% EndTime: 2019-03-09 08:58:43
% DurationCPUTime: 0.76s
% Computational Cost: add. (1794->203), mult. (4273->309), div. (0->0), fcn. (4866->12), ass. (0->103)
t169 = sin(pkin(12));
t172 = cos(pkin(12));
t176 = sin(qJ(5));
t220 = cos(qJ(5));
t154 = t176 * t169 - t220 * t172;
t175 = sin(qJ(6));
t178 = cos(qJ(6));
t189 = t178 * MDP(29) - t175 * MDP(30);
t186 = MDP(22) + t189;
t225 = t186 * t154;
t224 = 2 * MDP(15);
t223 = 0.2e1 * MDP(22);
t222 = 0.2e1 * MDP(29);
t221 = 0.2e1 * MDP(30);
t179 = cos(qJ(2));
t219 = pkin(1) * t179;
t218 = pkin(8) + qJ(3);
t170 = sin(pkin(11));
t161 = t170 * pkin(2) + qJ(4);
t217 = pkin(9) + t161;
t155 = t220 * t169 + t176 * t172;
t216 = t155 * t175;
t215 = t155 * t178;
t171 = sin(pkin(6));
t165 = t171 ^ 2;
t177 = sin(qJ(2));
t214 = t165 * t177;
t213 = t171 * t177;
t212 = t171 * t179;
t174 = cos(pkin(6));
t211 = t174 * MDP(8);
t159 = t174 * t219;
t141 = t174 * pkin(2) - t218 * t213 + t159;
t195 = t174 * t177 * pkin(1);
t145 = t218 * t212 + t195;
t173 = cos(pkin(11));
t132 = t170 * t141 + t173 * t145;
t127 = t174 * qJ(4) + t132;
t147 = t170 * t213 - t173 * t212;
t148 = (t170 * t179 + t173 * t177) * t171;
t156 = (-pkin(2) * t179 - pkin(1)) * t171;
t133 = t147 * pkin(3) - t148 * qJ(4) + t156;
t117 = t172 * t127 + t169 * t133;
t210 = t169 ^ 2 + t172 ^ 2;
t139 = t148 * t169 - t174 * t172;
t140 = t148 * t172 + t174 * t169;
t125 = -t176 * t139 + t220 * t140;
t119 = t175 * t125 - t147 * t178;
t209 = t119 * MDP(27);
t120 = t178 * t125 + t147 * t175;
t208 = t120 * MDP(24);
t207 = t120 * MDP(26);
t124 = t220 * t139 + t176 * t140;
t206 = t124 * MDP(28);
t205 = t125 * MDP(17);
t204 = t125 * MDP(18);
t203 = t147 * MDP(19);
t202 = t147 * MDP(20);
t201 = t154 * MDP(28);
t200 = t155 * MDP(23);
t163 = -t173 * pkin(2) - pkin(3);
t199 = t163 * MDP(16);
t198 = t169 * MDP(14);
t197 = t172 * MDP(13);
t196 = t178 * MDP(24);
t194 = t175 * t178 * MDP(25);
t193 = t217 * t169;
t116 = -t169 * t127 + t172 * t133;
t131 = t173 * t141 - t170 * t145;
t192 = t210 * MDP(16);
t191 = -t116 * t169 + t117 * t172;
t157 = -t172 * pkin(4) + t163;
t190 = MDP(26) * t178 - MDP(27) * t175;
t188 = MDP(29) * t175 + MDP(30) * t178;
t128 = -t174 * pkin(3) - t131;
t114 = t147 * pkin(4) - t140 * pkin(9) + t116;
t115 = -t139 * pkin(9) + t117;
t111 = t220 * t114 - t176 * t115;
t112 = t176 * t114 + t220 * t115;
t187 = -MDP(18) + t190;
t185 = (t177 * MDP(6) + t179 * MDP(7)) * t171;
t184 = -t197 + t198 + t200;
t121 = t139 * pkin(4) + t128;
t110 = t147 * pkin(10) + t112;
t113 = t124 * pkin(5) - t125 * pkin(10) + t121;
t107 = -t175 * t110 + t178 * t113;
t108 = t178 * t110 + t175 * t113;
t183 = t107 * MDP(29) - t108 * MDP(30) + t206 - t209;
t182 = t175 * MDP(26) + t178 * MDP(27) - t188 * pkin(10);
t181 = -MDP(20) + t182;
t168 = t178 ^ 2;
t167 = t175 ^ 2;
t152 = pkin(8) * t212 + t195;
t151 = -pkin(8) * t213 + t159;
t150 = t217 * t172;
t136 = t220 * t150 - t176 * t193;
t135 = t176 * t150 + t220 * t193;
t134 = t154 * pkin(5) - t155 * pkin(10) + t157;
t123 = t175 * t134 + t178 * t136;
t122 = t178 * t134 - t175 * t136;
t118 = t120 * t154;
t109 = -t147 * pkin(5) - t111;
t1 = [(t131 ^ 2 + t132 ^ 2 + t156 ^ 2) * MDP(12) + t147 ^ 2 * MDP(21) + (t116 ^ 2 + t117 ^ 2 + t128 ^ 2) * MDP(16) + MDP(1) + (MDP(4) * t177 + 0.2e1 * MDP(5) * t179) * t214 + (0.2e1 * t203 + t205) * t125 + (-0.2e1 * t119 * MDP(25) + t208) * t120 + (0.2e1 * t185 + t211) * t174 + (-0.2e1 * t202 - 0.2e1 * t204 + t206 + 0.2e1 * t207 - 0.2e1 * t209) * t124 + 0.2e1 * (t151 * t174 + t165 * t219) * MDP(9) + 0.2e1 * (-pkin(1) * t214 - t152 * t174) * MDP(10) + 0.2e1 * (-t117 * t147 + t128 * t140) * MDP(14) + 0.2e1 * (-t112 * t147 + t121 * t125) * MDP(23) + 0.2e1 * (-t131 * t148 - t132 * t147) * MDP(11) + (-t116 * t140 - t117 * t139) * t224 + 0.2e1 * (t116 * t147 + t128 * t139) * MDP(13) + (t111 * t147 + t121 * t124) * t223 + (t107 * t124 + t109 * t119) * t222 + (-t108 * t124 + t109 * t120) * t221; t211 + t151 * MDP(9) - t152 * MDP(10) + (-t128 * t172 + t163 * t139) * MDP(13) + (t128 * t169 + t163 * t140) * MDP(14) + t191 * MDP(15) + t128 * t199 + (t157 * t124 - t135 * t147) * MDP(22) + (t157 * t125 - t136 * t147) * MDP(23) + t118 * MDP(26) + (t135 * t119 + t122 * t124) * MDP(29) + (t135 * t120 - t123 * t124) * MDP(30) + t185 + ((-t139 * t172 + t140 * t169) * MDP(15) + t191 * MDP(16) + (-t169 * MDP(13) - t172 * MDP(14)) * t147) * t161 + (t121 * MDP(22) + t183 - t202 - t204) * t154 + ((-t147 * t170 - t148 * t173) * MDP(11) + (t131 * t173 + t132 * t170) * MDP(12)) * pkin(2) + (t121 * MDP(23) + (-t119 * t178 - t120 * t175) * MDP(25) + t120 * t196 + t205 + t203 + t188 * t109 + t187 * t124) * t155; 0.2e1 * t157 * t200 + MDP(8) + (t170 ^ 2 + t173 ^ 2) * MDP(12) * pkin(2) ^ 2 + (-0.2e1 * t197 + 0.2e1 * t198 + t199) * t163 + (t168 * MDP(24) + MDP(17) - 0.2e1 * t194) * t155 ^ 2 + (0.2e1 * t187 * t155 + t157 * t223 + t201) * t154 + (t122 * t154 + t135 * t216) * t222 + (-t123 * t154 + t135 * t215) * t221 + (t192 * t161 + t210 * t224) * t161; t156 * MDP(12) + (-t169 * t139 - t172 * t140) * MDP(15) + (t116 * t172 + t117 * t169) * MDP(16) + (t154 * t119 - t124 * t216) * MDP(29) + (-t124 * t215 + t118) * MDP(30) + (-t154 * MDP(22) - t184) * t147; 0; MDP(12) + t192; t139 * MDP(13) + t140 * MDP(14) + t128 * MDP(16) + t125 * MDP(23) + t186 * t124; t184 + t199 + t225; 0; MDP(16); t125 * MDP(19) + t147 * MDP(21) + t111 * MDP(22) - t112 * MDP(23) + t175 * t208 + (-t175 * t119 + t120 * t178) * MDP(25) + (-pkin(5) * t119 - t109 * t178) * MDP(29) + (-pkin(5) * t120 + t109 * t175) * MDP(30) + t181 * t124; -t136 * MDP(23) - t186 * t135 + t181 * t154 + (MDP(19) + t175 * t196 + (-t167 + t168) * MDP(25) - t188 * pkin(5)) * t155; -t200 - t225; 0; t167 * MDP(24) + 0.2e1 * pkin(5) * t189 + MDP(21) + 0.2e1 * t194; t183 + t207; t122 * MDP(29) - t123 * MDP(30) + t190 * t155 + t201; -t188 * t155; t189; t182; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
