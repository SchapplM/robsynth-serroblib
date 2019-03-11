% Calculate joint inertia matrix for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR12_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:21:50
% EndTime: 2019-03-09 11:21:54
% DurationCPUTime: 1.20s
% Computational Cost: add. (1344->216), mult. (2835->319), div. (0->0), fcn. (3057->10), ass. (0->89)
t171 = sin(pkin(6));
t175 = sin(qJ(2));
t209 = t171 * t175;
t154 = pkin(8) * t209;
t173 = cos(pkin(6));
t178 = cos(qJ(2));
t213 = pkin(1) * t178;
t197 = -pkin(2) - t213;
t132 = pkin(3) * t209 + t154 + (-pkin(9) + t197) * t173;
t179 = -pkin(2) - pkin(9);
t194 = -qJ(3) * t175 - pkin(1);
t137 = (t179 * t178 + t194) * t171;
t177 = cos(qJ(4));
t214 = sin(qJ(4));
t123 = t177 * t132 - t214 * t137;
t124 = t214 * t132 + t177 * t137;
t208 = t171 * t178;
t142 = t173 * t214 + t177 * t208;
t143 = t173 * t177 - t214 * t208;
t223 = t143 * MDP(17) - t142 * MDP(18) + t123 * MDP(20) - t124 * MDP(21);
t170 = sin(pkin(11));
t172 = cos(pkin(11));
t149 = -t170 * t214 + t172 * t177;
t174 = sin(qJ(6));
t176 = cos(qJ(6));
t221 = (MDP(26) * t176 - MDP(27) * t174) * t149;
t220 = -qJ(5) + t179;
t198 = t176 * MDP(30);
t199 = t174 * MDP(29);
t185 = t198 + t199;
t219 = t185 * (t170 * pkin(4) + pkin(10)) - t174 * MDP(26) - t176 * MDP(27);
t218 = -2 * MDP(16);
t217 = 2 * MDP(22);
t216 = 0.2e1 * MDP(29);
t215 = 0.2e1 * MDP(30);
t212 = pkin(2) * MDP(14);
t152 = t220 * t214;
t193 = t220 * t177;
t133 = t170 * t152 - t172 * t193;
t211 = t133 * t149;
t128 = t172 * t142 + t170 * t143;
t150 = -t170 * t177 - t172 * t214;
t210 = t150 * t128;
t207 = t173 * t178;
t163 = t214 * pkin(4) + qJ(3);
t118 = pkin(4) * t209 - t143 * qJ(5) + t123;
t120 = -t142 * qJ(5) + t124;
t115 = t170 * t118 + t172 * t120;
t145 = t173 * t175 * pkin(1) + pkin(8) * t208;
t129 = -t170 * t142 + t172 * t143;
t125 = t174 * t129 - t176 * t209;
t206 = t125 * MDP(25);
t205 = t125 * MDP(27);
t126 = t176 * t129 + t174 * t209;
t204 = t126 * MDP(24);
t203 = t128 * MDP(26);
t202 = t128 * MDP(28);
t201 = t143 * MDP(15);
t200 = t174 * MDP(24);
t164 = t173 * qJ(3);
t138 = -t164 - t145;
t196 = t174 * t176 * MDP(25);
t195 = t214 * MDP(21);
t136 = pkin(3) * t208 - t138;
t114 = t172 * t118 - t170 * t120;
t192 = t149 * t172 - t150 * t170;
t191 = -t145 * MDP(10) + (pkin(1) * t207 - t154) * MDP(9);
t186 = t176 * MDP(29) - t174 * MDP(30);
t127 = t142 * pkin(4) + t136;
t184 = t177 * MDP(20) - t195;
t113 = pkin(10) * t209 + t115;
t116 = t128 * pkin(5) - t129 * pkin(10) + t127;
t110 = -t174 * t113 + t176 * t116;
t111 = t176 * t113 + t174 * t116;
t182 = t126 * MDP(26) + t110 * MDP(29) - t111 * MDP(30) + t202 - t205;
t181 = -t179 * t195 - t214 * MDP(18) + (t179 * MDP(20) + MDP(17)) * t177;
t169 = t176 ^ 2;
t167 = t174 ^ 2;
t162 = -t172 * pkin(4) - pkin(5);
t148 = t150 ^ 2;
t147 = t149 ^ 2;
t140 = t197 * t173 + t154;
t139 = (-pkin(2) * t178 + t194) * t171;
t135 = t172 * t152 + t170 * t193;
t130 = -t150 * pkin(5) - t149 * pkin(10) + t163;
t122 = t174 * t130 + t176 * t135;
t121 = t176 * t130 - t174 * t135;
t112 = -pkin(5) * t209 - t114;
t1 = [(t114 ^ 2 + t115 ^ 2 + t127 ^ 2) * MDP(23) + MDP(1) + (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) * MDP(14) + t173 ^ 2 * MDP(8) + (t142 * t218 + t201) * t143 + (t202 - 0.2e1 * t205) * t128 + (0.2e1 * t203 + t204 - 0.2e1 * t206) * t126 + (t110 * t128 + t112 * t125) * t216 + (-t111 * t128 + t112 * t126) * t215 + (-t114 * t129 - t115 * t128) * t217 + 0.2e1 * (t140 * MDP(12) - t138 * MDP(13) + t191) * t173 + 0.2e1 * (t142 * MDP(20) + t143 * MDP(21)) * t136 + ((0.2e1 * t178 * MDP(5) + (MDP(19) + MDP(4)) * t175) * t175 + 0.2e1 * (-MDP(10) * t175 + MDP(9) * t178) * pkin(1)) * t171 ^ 2 + 0.2e1 * (MDP(7) * t207 + (-t138 * MDP(11) + t139 * MDP(12)) * t178 + (t140 * MDP(11) - t139 * MDP(13) + MDP(6) * t173 + t223) * t175) * t171; t154 * MDP(12) + (0.2e1 * t164 + t145) * MDP(13) + (-t140 * pkin(2) - t138 * qJ(3)) * MDP(14) + (-t177 * t142 - t143 * t214) * MDP(16) + (qJ(3) * t142 + t136 * t214) * MDP(20) + (qJ(3) * t143 + t136 * t177) * MDP(21) + (-t135 * t128 + t133 * t129) * MDP(22) + (-t114 * t133 + t115 * t135 + t127 * t163) * MDP(23) + (t121 * t128 + t133 * t125) * MDP(29) + (-t122 * t128 + t133 * t126) * MDP(30) + t177 * t201 + (MDP(8) + (-0.2e1 * pkin(2) - t213) * MDP(12)) * t173 + (t115 * MDP(22) - t182) * t150 + (-t114 * MDP(22) + (-t126 * MDP(25) - t128 * MDP(27) + t112 * MDP(29)) * t174 + (t112 * MDP(30) + t203 + t204 - t206) * t176) * t149 + ((qJ(3) * MDP(11) + MDP(7)) * t178 + (-pkin(2) * MDP(11) + MDP(6) + t181) * t175) * t171 + t191; MDP(8) + (t133 ^ 2 + t135 ^ 2 + t163 ^ 2) * MDP(23) + t148 * MDP(28) + (MDP(15) * t177 + t214 * t218) * t177 - 0.2e1 * t150 * t221 + (t169 * MDP(24) - 0.2e1 * t196) * t147 + (-0.2e1 * MDP(12) + t212) * pkin(2) + (MDP(14) * qJ(3) + 0.2e1 * t214 * MDP(20) + 0.2e1 * t177 * MDP(21) + 0.2e1 * MDP(13)) * qJ(3) + (t135 * t150 + t211) * t217 + (-t121 * t150 + t174 * t211) * t216 + (t122 * t150 + t176 * t211) * t215; t173 * MDP(12) + t140 * MDP(14) + (-t149 * t129 + t210) * MDP(22) + (t114 * t149 - t115 * t150) * MDP(23) + (-t149 * t125 + t174 * t210) * MDP(29) + (-t149 * t126 + t176 * t210) * MDP(30) + (MDP(11) + t184) * t209; -t148 * t199 - t212 + MDP(12) + (-t135 * MDP(23) + (-MDP(22) - t198) * t150) * t150 + (-t133 * MDP(23) + (-MDP(22) - t185) * t149) * t149; MDP(14) + (t148 + t147) * MDP(23); MDP(19) * t209 + t126 * t200 + (-t174 * t125 + t126 * t176) * MDP(25) + (-t112 * t176 + t162 * t125) * MDP(29) + (t112 * t174 + t162 * t126) * MDP(30) - t219 * t128 + ((-t128 * t170 - t129 * t172) * MDP(22) + (t114 * t172 + t115 * t170) * MDP(23)) * pkin(4) + t223; -t186 * t133 + t219 * t150 + (t176 * t200 + (-t167 + t169) * MDP(25) + t185 * t162) * t149 + (-t192 * MDP(22) + (-t133 * t172 + t135 * t170) * MDP(23)) * pkin(4) + t181; MDP(23) * pkin(4) * t192 + t149 * t186 + t184; 0.2e1 * t196 + t167 * MDP(24) + MDP(19) + (t170 ^ 2 + t172 ^ 2) * MDP(23) * pkin(4) ^ 2 - 0.2e1 * t186 * t162; t127 * MDP(23) + t128 * t186; t163 * MDP(23) - t150 * t186; 0; 0; MDP(23); t182; -t150 * MDP(28) + t121 * MDP(29) - t122 * MDP(30) + t221; t185 * t150; -t219; t186; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
