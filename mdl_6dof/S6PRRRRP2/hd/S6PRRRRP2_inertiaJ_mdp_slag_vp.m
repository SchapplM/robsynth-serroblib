% Calculate joint inertia matrix for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRRP2_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:04:22
% EndTime: 2019-03-09 00:04:24
% DurationCPUTime: 0.72s
% Computational Cost: add. (892->174), mult. (1758->242), div. (0->0), fcn. (1920->10), ass. (0->94)
t159 = sin(qJ(4));
t145 = t159 * pkin(3) + pkin(10);
t158 = sin(qJ(5));
t154 = t158 ^ 2;
t162 = cos(qJ(5));
t155 = t162 ^ 2;
t205 = t154 + t155;
t207 = t205 * t145;
t225 = t158 * MDP(21) + t162 * MDP(22);
t196 = t158 * MDP(25);
t224 = t196 - MDP(17);
t160 = sin(qJ(3));
t217 = cos(qJ(4));
t218 = cos(qJ(3));
t135 = t159 * t218 + t217 * t160;
t223 = 0.2e1 * t135;
t222 = pkin(9) + pkin(8);
t221 = MDP(24) + MDP(26);
t147 = -t218 * pkin(3) - pkin(2);
t220 = 0.2e1 * t147;
t219 = 2 * MDP(27);
t134 = t159 * t160 - t217 * t218;
t216 = pkin(10) * t134;
t215 = t134 * pkin(5);
t193 = t217 * pkin(3);
t146 = -t193 - pkin(4);
t214 = pkin(4) - t146;
t213 = t134 * qJ(6);
t212 = t134 * t145;
t156 = sin(pkin(6));
t161 = sin(qJ(2));
t211 = t156 * t161;
t163 = cos(qJ(2));
t210 = t156 * t163;
t209 = t158 * t162;
t118 = t134 * pkin(4) - t135 * pkin(10) + t147;
t138 = t222 * t160;
t139 = t222 * t218;
t124 = -t159 * t138 + t217 * t139;
t105 = t158 * t118 + t162 * t124;
t181 = -t162 * pkin(5) - t158 * qJ(6);
t137 = -pkin(4) + t181;
t132 = -t193 + t137;
t208 = -t132 - t137;
t206 = t205 * pkin(10);
t204 = MDP(18) * t135;
t203 = MDP(29) * t158;
t102 = t105 + t213;
t202 = t102 * MDP(29);
t186 = -t162 * t118 + t158 * t124;
t103 = t186 - t215;
t201 = t103 * MDP(29);
t123 = t217 * t138 + t159 * t139;
t180 = pkin(5) * t158 - t162 * qJ(6);
t106 = t180 * t135 + t123;
t200 = t106 * MDP(28);
t199 = t132 * MDP(29);
t198 = t134 * MDP(23);
t197 = t137 * MDP(29);
t195 = t158 * MDP(28);
t194 = MDP(25) - MDP(28);
t192 = -t180 * MDP(27) + t225;
t191 = MDP(20) * t209;
t190 = t154 * MDP(19) + MDP(16) + 0.2e1 * t191;
t189 = MDP(10) * t218;
t188 = -MDP(29) * pkin(5) - MDP(26);
t185 = t205 * MDP(29);
t157 = cos(pkin(6));
t128 = t157 * t218 - t160 * t211;
t129 = t157 * t160 + t218 * t211;
t113 = -t217 * t128 + t159 * t129;
t114 = t159 * t128 + t217 * t129;
t108 = t158 * t114 + t162 * t210;
t109 = t162 * t114 - t158 * t210;
t177 = t108 * t158 + t109 * t162;
t184 = -t114 * MDP(18) + t177 * MDP(27) + t224 * t113;
t183 = -pkin(4) * t135 - t216;
t182 = -MDP(24) + t188;
t179 = -t135 * t137 + t216;
t178 = MDP(29) * qJ(6) - t194;
t176 = t132 * t135 - t212;
t175 = t135 * t146 - t212;
t174 = t162 * MDP(21) - t158 * MDP(22);
t173 = -t186 * MDP(24) - t105 * MDP(25);
t172 = t162 * MDP(24) - t196;
t171 = -t123 * MDP(24) - t106 * MDP(26);
t170 = -0.2e1 * t162 * MDP(26) - 0.2e1 * t195;
t169 = t177 * MDP(29);
t168 = -t221 * t162 - t195;
t167 = (t217 * MDP(17) - t159 * MDP(18)) * pkin(3);
t166 = (t102 * t162 + t103 * t158) * MDP(27) - t124 * MDP(18) + t224 * t123 + ((-t154 + t155) * MDP(20) + MDP(19) * t209 + MDP(14)) * t135 + (-MDP(15) + t225) * t134;
t165 = t182 * t158 + t178 * t162;
t149 = t158 * MDP(27);
t1 = [MDP(1) + (t108 ^ 2 + t109 ^ 2 + t113 ^ 2) * MDP(29); (t108 * t103 + t113 * t106) * MDP(29) + (-t194 * t134 + t202) * t109 + ((t108 * t162 - t109 * t158) * MDP(27) + t194 * t162 * t113) * t135 + (-t161 * MDP(4) + (-MDP(11) * t160 - MDP(17) * t134 + MDP(3) + t189 - t204) * t163) * t156 + t221 * (t113 * t158 * t135 - t108 * t134); MDP(2) + 0.2e1 * pkin(2) * t189 + t204 * t220 + (t102 ^ 2 + t103 ^ 2 + t106 ^ 2) * MDP(29) + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t160 + 0.2e1 * t218 * MDP(6)) * t160 + (t155 * MDP(19) + MDP(12) - 0.2e1 * t191) * t135 ^ 2 + (MDP(17) * t220 + t198 + (-MDP(13) + t174) * t223) * t134 + 0.2e1 * (-t103 * MDP(26) + t102 * MDP(28) + t173) * t134 + ((-t102 * t158 + t103 * t162) * MDP(27) + (t158 * MDP(24) + t162 * MDP(25)) * t123 + (t158 * MDP(26) - t162 * MDP(28)) * t106) * t223; t128 * MDP(10) - t129 * MDP(11) + t145 * t169 + (t168 + t199) * t113 + t184; (t175 * MDP(24) + t176 * MDP(26) + t145 * t201 - t200) * t158 + t106 * t199 + (t175 * MDP(25) - t176 * MDP(28) + t145 * t202 + t171) * t162 + t160 * MDP(7) + t218 * MDP(8) + (-t160 * MDP(10) - t218 * MDP(11)) * pkin(8) + t166; MDP(9) + t207 * t219 + t145 ^ 2 * t185 + (t170 + t199) * t132 + t190 - 0.2e1 * t172 * t146 + 0.2e1 * t167; pkin(10) * t169 + (t168 + t197) * t113 + t184; t106 * t197 + (t183 * MDP(25) + t179 * MDP(28) + pkin(10) * t202 + t171) * t162 + (t183 * MDP(24) - t179 * MDP(26) + pkin(10) * t201 - t200) * t158 + t166; (t206 + t207) * MDP(27) + (pkin(10) * t207 + t132 * t137) * MDP(29) + t167 + (t214 * MDP(24) + t208 * MDP(26)) * t162 + (-t214 * MDP(25) + t208 * MDP(28)) * t158 + t190; t206 * t219 + pkin(10) ^ 2 * t185 + (t170 + t197) * t137 + 0.2e1 * t172 * pkin(4) + t190; t182 * t108 + t178 * t109; t198 + (-t186 + 0.2e1 * t215) * MDP(26) + (t105 + 0.2e1 * t213) * MDP(28) + (-t103 * pkin(5) + t102 * qJ(6)) * MDP(29) + (t181 * MDP(27) + t174) * t135 + t173; t165 * t145 + t192; t165 * pkin(10) + t192; MDP(23) + 0.2e1 * pkin(5) * MDP(26) + 0.2e1 * qJ(6) * MDP(28) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(29); t108 * MDP(29); t162 * t135 * MDP(27) - t134 * MDP(26) + t201; t145 * t203 + t149; pkin(10) * t203 + t149; t188; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
