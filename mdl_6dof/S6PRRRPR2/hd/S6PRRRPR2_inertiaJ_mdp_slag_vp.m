% Calculate joint inertia matrix for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR2_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:09:00
% EndTime: 2019-03-08 23:09:02
% DurationCPUTime: 0.72s
% Computational Cost: add. (912->171), mult. (1859->247), div. (0->0), fcn. (2122->12), ass. (0->91)
t186 = sin(qJ(4));
t170 = t186 * pkin(3) + qJ(5);
t181 = sin(pkin(12));
t183 = cos(pkin(12));
t217 = t181 ^ 2 + t183 ^ 2;
t219 = t217 * t170;
t234 = t183 * MDP(19) - t181 * MDP(20);
t185 = sin(qJ(6));
t189 = cos(qJ(6));
t157 = t185 * t181 - t189 * t183;
t158 = t189 * t181 + t185 * t183;
t235 = t157 * MDP(28) + t158 * MDP(29);
t198 = -t234 + t235;
t221 = t158 * MDP(25) - t157 * MDP(26);
t190 = cos(qJ(3));
t174 = -t190 * pkin(3) - pkin(2);
t233 = 0.2e1 * t174;
t232 = 2 * MDP(21);
t231 = -2 * MDP(24);
t230 = -pkin(9) - pkin(8);
t229 = cos(qJ(4));
t228 = t183 * pkin(5);
t178 = t183 * pkin(10);
t226 = pkin(4) * MDP(22);
t187 = sin(qJ(3));
t166 = t230 * t187;
t167 = t230 * t190;
t144 = -t229 * t166 - t186 * t167;
t225 = t144 * t183;
t160 = t186 * t190 + t229 * t187;
t224 = t160 * t181;
t182 = sin(pkin(6));
t188 = sin(qJ(2));
t223 = t182 * t188;
t191 = cos(qJ(2));
t222 = t182 * t191;
t159 = t186 * t187 - t229 * t190;
t134 = t159 * pkin(4) - t160 * qJ(5) + t174;
t145 = t186 * t166 - t229 * t167;
t114 = t181 * t134 + t183 * t145;
t218 = t217 * qJ(5);
t216 = MDP(10) * t190;
t129 = t157 * t160;
t215 = MDP(23) * t129;
t128 = t158 * t160;
t214 = t128 * MDP(26);
t213 = t129 * MDP(25);
t212 = t159 * MDP(27);
t173 = -t229 * pkin(3) - pkin(4);
t211 = t173 * MDP(22);
t209 = MDP(16) + (MDP(23) * t158 + t157 * t231) * t158;
t113 = t183 * t134 - t181 * t145;
t208 = t217 * MDP(22);
t207 = -pkin(4) * t160 - qJ(5) * t159;
t206 = -t113 * t181 + t114 * t183;
t184 = cos(pkin(6));
t146 = t184 * t190 - t187 * t223;
t147 = t184 * t187 + t190 * t223;
t127 = t186 * t146 + t229 * t147;
t120 = -t181 * t127 - t183 * t222;
t121 = t183 * t127 - t181 * t222;
t205 = -t120 * t181 + t121 * t183;
t204 = -t159 * t170 + t160 * t173;
t203 = 0.2e1 * t234;
t202 = t181 * MDP(19) + t183 * MDP(20);
t111 = t159 * pkin(5) - t160 * t178 + t113;
t112 = -pkin(10) * t224 + t114;
t201 = (t189 * t111 - t185 * t112) * MDP(28) - (t185 * t111 + t189 * t112) * MDP(29);
t200 = (t189 * t120 - t185 * t121) * MDP(28) - (t185 * t120 + t189 * t121) * MDP(29);
t199 = t128 * MDP(28) - t129 * MDP(29);
t197 = 0.2e1 * t235;
t196 = (t229 * MDP(17) - t186 * MDP(18)) * pkin(3);
t195 = t144 * MDP(22) + t199;
t194 = t206 * MDP(21) + (-t158 * t128 + t129 * t157) * MDP(24) - t158 * t215 - t144 * MDP(17) - t145 * MDP(18) + t160 * MDP(14) + (-MDP(15) + t221) * t159;
t126 = -t229 * t146 + t186 * t147;
t193 = -t127 * MDP(18) + t205 * MDP(21) + (-MDP(17) + t198) * t126;
t171 = -pkin(4) - t228;
t163 = t183 * qJ(5) + t178;
t162 = (-pkin(10) - qJ(5)) * t181;
t161 = t173 - t228;
t154 = t183 * t170 + t178;
t153 = (-pkin(10) - t170) * t181;
t141 = t185 * t162 + t189 * t163;
t140 = t189 * t162 - t185 * t163;
t137 = t144 * t181;
t133 = t185 * t153 + t189 * t154;
t132 = t189 * t153 - t185 * t154;
t123 = pkin(5) * t224 + t144;
t116 = t123 * t158;
t115 = t123 * t157;
t1 = [MDP(1) + (t120 ^ 2 + t121 ^ 2 + t126 ^ 2) * MDP(22); (t120 * t113 + t121 * t114) * MDP(22) + t195 * t126 + (t120 * MDP(19) - t121 * MDP(20) + t200) * t159 + ((-t120 * t183 - t121 * t181) * MDP(21) + t202 * t126) * t160 + (-t188 * MDP(4) + (-MDP(11) * t187 - MDP(17) * t159 - MDP(18) * t160 + MDP(3) + t216) * t191) * t182; MDP(2) + 0.2e1 * pkin(2) * t216 + (t113 ^ 2 + t114 ^ 2 + t144 ^ 2) * MDP(22) + (MDP(12) * t160 + MDP(18) * t233) * t160 - (t128 * t231 - t215) * t129 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t187 + 0.2e1 * t190 * MDP(6)) * t187 + (-0.2e1 * t160 * MDP(13) + MDP(17) * t233 + t212 - 0.2e1 * t213 - 0.2e1 * t214) * t159 + 0.2e1 * t199 * t123 + 0.2e1 * (t113 * MDP(19) - t114 * MDP(20) + t201) * t159 + 0.2e1 * ((-t113 * t183 - t114 * t181) * MDP(21) + t202 * t144) * t160; t146 * MDP(10) - t147 * MDP(11) + (t126 * t173 + t205 * t170) * MDP(22) + t193; (t204 * t181 - t225) * MDP(19) + (t204 * t183 + t137) * MDP(20) + (t144 * t173 + t206 * t170) * MDP(22) + (-t187 * MDP(10) - t190 * MDP(11)) * pkin(8) + t194 + (t161 * t128 + t132 * t159 + t115) * MDP(28) + (-t161 * t129 - t133 * t159 + t116) * MDP(29) + t187 * MDP(7) + t190 * MDP(8); MDP(9) + t219 * t232 + t170 ^ 2 * t208 + (-t203 + t211) * t173 + t161 * t197 + 0.2e1 * t196 + t209; (-t126 * pkin(4) + t205 * qJ(5)) * MDP(22) + t193; (t207 * t181 - t225) * MDP(19) + (t207 * t183 + t137) * MDP(20) + (-t144 * pkin(4) + t206 * qJ(5)) * MDP(22) + (t171 * t128 + t140 * t159 + t115) * MDP(28) + (-t171 * t129 - t141 * t159 + t116) * MDP(29) + t194; (t218 + t219) * MDP(21) + (-t173 * pkin(4) + qJ(5) * t219) * MDP(22) + t196 + t209 + t234 * (pkin(4) - t173) + t235 * (t161 + t171); t218 * t232 + qJ(5) ^ 2 * t208 + t171 * t197 + (t203 + t226) * pkin(4) + t209; t126 * MDP(22); t202 * t160 + t195; t198 + t211; t198 - t226; MDP(22); t200; t201 + t212 - t213 - t214; t132 * MDP(28) - t133 * MDP(29) + t221; t140 * MDP(28) - t141 * MDP(29) + t221; 0; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
