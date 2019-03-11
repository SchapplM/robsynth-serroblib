% Calculate joint inertia matrix for
% S6RRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP2_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:45:15
% EndTime: 2019-03-09 11:45:17
% DurationCPUTime: 0.70s
% Computational Cost: add. (1517->170), mult. (2744->238), div. (0->0), fcn. (3146->8), ass. (0->85)
t157 = cos(pkin(10));
t146 = t157 * pkin(2) + pkin(3);
t159 = sin(qJ(4));
t156 = sin(pkin(10));
t204 = t156 * pkin(2);
t207 = cos(qJ(4));
t134 = -t159 * t146 - t207 * t204;
t132 = pkin(9) - t134;
t158 = sin(qJ(5));
t154 = t158 ^ 2;
t161 = cos(qJ(5));
t155 = t161 ^ 2;
t195 = t154 + t155;
t197 = t195 * t132;
t212 = t158 * MDP(22) + t161 * MDP(23);
t160 = sin(qJ(2));
t162 = cos(qJ(2));
t136 = t156 * t162 + t157 * t160;
t171 = t156 * t160 - t157 * t162;
t123 = t207 * t136 - t159 * t171;
t211 = 0.2e1 * t123;
t147 = -t162 * pkin(2) - pkin(1);
t129 = t171 * pkin(3) + t147;
t210 = 0.2e1 * t129;
t209 = 0.2e1 * t162;
t208 = 2 * MDP(28);
t122 = t159 * t136 + t207 * t171;
t206 = pkin(9) * t122;
t205 = t122 * pkin(5);
t133 = t207 * t146 - t159 * t204;
t131 = -pkin(4) - t133;
t203 = pkin(4) - t131;
t202 = -qJ(3) - pkin(7);
t201 = t122 * qJ(6);
t200 = t122 * t132;
t199 = t158 * t161;
t110 = t122 * pkin(4) - t123 * pkin(9) + t129;
t142 = t202 * t160;
t143 = t202 * t162;
t125 = t157 * t142 + t156 * t143;
t114 = -t136 * pkin(8) + t125;
t126 = t156 * t142 - t157 * t143;
t115 = -t171 * pkin(8) + t126;
t112 = t159 * t114 + t207 * t115;
t102 = t158 * t110 + t161 * t112;
t179 = -t161 * pkin(5) - t158 * qJ(6);
t141 = -pkin(4) + t179;
t124 = t141 - t133;
t198 = -t124 - t141;
t196 = t195 * pkin(9);
t194 = MDP(30) * t141;
t193 = MDP(30) * t158;
t183 = -t161 * t110 + t158 * t112;
t100 = t183 - t205;
t192 = t100 * MDP(30);
t191 = t122 * MDP(24);
t190 = t133 * MDP(18);
t189 = t134 * MDP(19);
t188 = -MDP(26) + MDP(29);
t178 = pkin(5) * t158 - t161 * qJ(6);
t187 = -t178 * MDP(28) + t212;
t186 = MDP(21) * t199;
t185 = t154 * MDP(20) + MDP(17) + 0.2e1 * t186;
t184 = -MDP(30) * pkin(5) - MDP(27);
t182 = t195 * MDP(30);
t181 = -pkin(4) * t123 - t206;
t180 = MDP(25) - t184;
t177 = -t123 * t141 + t206;
t99 = t102 + t201;
t176 = t100 * t161 - t99 * t158;
t175 = t100 * t158 + t99 * t161;
t174 = MDP(30) * qJ(6) + t188;
t173 = -t123 * t124 + t200;
t172 = t123 * t131 - t200;
t170 = t161 * MDP(22) - t158 * MDP(23);
t169 = -t183 * MDP(25) - t102 * MDP(26);
t168 = t161 * MDP(25) - t158 * MDP(26);
t167 = -0.2e1 * t161 * MDP(27) - 0.2e1 * t158 * MDP(29);
t111 = -t207 * t114 + t159 * t115;
t166 = -t111 * MDP(18) - t112 * MDP(19) + t175 * MDP(28) + ((-t154 + t155) * MDP(21) + MDP(20) * t199 + MDP(15)) * t123 + (-MDP(16) + t212) * t122;
t165 = -t180 * t158 + t174 * t161;
t149 = t158 * MDP(28);
t107 = t111 * t158;
t103 = t178 * t123 + t111;
t1 = [MDP(1) + pkin(1) * MDP(9) * t209 + (t125 ^ 2 + t126 ^ 2 + t147 ^ 2) * MDP(12) + (t100 ^ 2 + t103 ^ 2 + t99 ^ 2) * MDP(30) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t160 + MDP(5) * t209) * t160 + (MDP(18) * t210 + t191 + (-MDP(14) + t170) * t211) * t122 + 0.2e1 * (-t125 * t136 - t126 * t171) * MDP(11) + 0.2e1 * (-t100 * MDP(27) + t99 * MDP(29) + t169) * t122 + (t176 * MDP(28) + (t158 * MDP(25) + t161 * MDP(26)) * t111 + (t158 * MDP(27) - t161 * MDP(29)) * t103) * t211 + (MDP(19) * t210 + (t155 * MDP(20) + MDP(13) - 0.2e1 * t186) * t123) * t123; t162 * MDP(7) + t160 * MDP(6) + t166 + (t172 * t161 + t107) * MDP(26) + (t103 * t124 + t175 * t132) * MDP(30) + (-t103 * t158 + t173 * t161) * MDP(29) + (-t103 * t161 - t173 * t158) * MDP(27) + (-t111 * t161 + t172 * t158) * MDP(25) + (-t162 * MDP(10) - t160 * MDP(9)) * pkin(7) + ((t125 * t157 + t126 * t156) * MDP(12) + (-t157 * t136 - t156 * t171) * MDP(11)) * pkin(2); MDP(8) + (t156 ^ 2 + t157 ^ 2) * MDP(12) * pkin(2) ^ 2 - 0.2e1 * t168 * t131 + t132 ^ 2 * t182 + (t124 * MDP(30) + t167) * t124 + 0.2e1 * t190 + 0.2e1 * t189 + t197 * t208 + t185; t147 * MDP(12) - t176 * MDP(30) + (-t195 * MDP(28) + MDP(19)) * t123 + (MDP(18) + (MDP(25) + MDP(27)) * t161 + t188 * t158) * t122; 0; MDP(12) + t182; t103 * t194 + t107 * MDP(26) + (pkin(9) * t99 * MDP(30) - t111 * MDP(25) + t181 * MDP(26) - t103 * MDP(27) + t177 * MDP(29)) * t161 + (t181 * MDP(25) - t177 * MDP(27) - t103 * MDP(29) + pkin(9) * t192) * t158 + t166; t190 + t189 + (t196 + t197) * MDP(28) + (pkin(9) * t197 + t124 * t141) * MDP(30) + (t203 * MDP(25) + t198 * MDP(27)) * t161 + (-t203 * MDP(26) + t198 * MDP(29)) * t158 + t185; 0; t196 * t208 + pkin(9) ^ 2 * t182 + (t167 + t194) * t141 + 0.2e1 * t168 * pkin(4) + t185; t191 + (-t183 + 0.2e1 * t205) * MDP(27) + (t102 + 0.2e1 * t201) * MDP(29) + (-t100 * pkin(5) + t99 * qJ(6)) * MDP(30) + (t179 * MDP(28) + t170) * t123 + t169; t165 * t132 + t187; t174 * t158 + t180 * t161; t165 * pkin(9) + t187; MDP(24) + 0.2e1 * pkin(5) * MDP(27) + 0.2e1 * qJ(6) * MDP(29) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(30); t161 * t123 * MDP(28) - t122 * MDP(27) + t192; t132 * t193 + t149; -t161 * MDP(30); pkin(9) * t193 + t149; t184; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
