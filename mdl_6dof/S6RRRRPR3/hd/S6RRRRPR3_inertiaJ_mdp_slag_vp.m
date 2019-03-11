% Calculate joint inertia matrix for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR3_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:04:02
% EndTime: 2019-03-09 22:04:03
% DurationCPUTime: 0.57s
% Computational Cost: add. (1088->180), mult. (1940->226), div. (0->0), fcn. (2120->8), ass. (0->92)
t216 = MDP(23) - MDP(26);
t189 = -MDP(24) + MDP(27);
t159 = sin(qJ(4));
t215 = t189 * t159;
t158 = sin(qJ(6));
t162 = cos(qJ(6));
t207 = 2 * MDP(35);
t208 = 2 * MDP(34);
t172 = t158 * t208 + t162 * t207 + (2 * MDP(27));
t164 = cos(qJ(3));
t150 = pkin(2) * t164 + pkin(3);
t163 = cos(qJ(4));
t160 = sin(qJ(3));
t204 = pkin(2) * t160;
t196 = -t163 * t150 + t159 * t204;
t129 = -pkin(4) + t196;
t214 = MDP(28) * t129;
t165 = cos(qJ(2));
t151 = -t165 * pkin(2) - pkin(1);
t213 = 0.2e1 * t151;
t212 = 0.2e1 * t165;
t211 = -2 * MDP(26);
t210 = 2 * MDP(26);
t206 = pkin(4) + pkin(10);
t205 = pkin(7) + pkin(8);
t203 = pkin(3) * t163;
t201 = MDP(28) * pkin(4);
t161 = sin(qJ(2));
t135 = t160 * t161 - t164 * t165;
t136 = t160 * t165 + t161 * t164;
t119 = t163 * t135 + t136 * t159;
t200 = qJ(5) * t119;
t140 = t205 * t161;
t141 = t205 * t165;
t185 = -t164 * t140 - t141 * t160;
t112 = -pkin(9) * t136 + t185;
t178 = t140 * t160 - t141 * t164;
t113 = -pkin(9) * t135 - t178;
t108 = t112 * t159 + t113 * t163;
t101 = -pkin(5) * t119 + t108;
t97 = t101 * t158;
t98 = t101 * t162;
t174 = t150 * t159 + t163 * t204;
t127 = qJ(5) + t174;
t199 = t119 * t127;
t146 = pkin(3) * t159 + qJ(5);
t198 = t119 * t146;
t197 = t158 * t162;
t149 = -pkin(4) - t203;
t195 = MDP(28) * t149;
t194 = MDP(34) * t162;
t120 = -t135 * t159 + t136 * t163;
t193 = t120 * MDP(32);
t192 = t196 * MDP(23);
t191 = t174 * MDP(24);
t153 = t162 * MDP(31);
t190 = t164 * MDP(16);
t188 = 0.2e1 * t120;
t187 = MDP(30) * t197;
t157 = t162 ^ 2;
t186 = t157 * MDP(29) + MDP(22) - 0.2e1 * t187;
t184 = -t158 * MDP(32) + t153;
t183 = MDP(15) + t186;
t182 = -t196 + t203;
t181 = t120 * t206 + t200;
t107 = -t112 * t163 + t113 * t159;
t126 = -pkin(10) + t129;
t180 = -t120 * t126 + t199;
t145 = -pkin(10) + t149;
t179 = -t120 * t145 + t198;
t125 = t135 * pkin(3) + t151;
t177 = MDP(31) * t158 + MDP(32) * t162;
t176 = -MDP(35) * t158 + t194;
t167 = 0.2e1 * qJ(5);
t173 = t167 * MDP(27) + t186;
t171 = -t120 * qJ(5) + t125;
t156 = t158 ^ 2;
t170 = t98 * MDP(35) + (t153 + MDP(20)) * t120 + t189 * t108 - t216 * t107 + ((-t156 + t157) * MDP(30) + MDP(29) * t197 - MDP(21)) * t119;
t169 = t136 * MDP(13) - t135 * MDP(14) + t185 * MDP(16) + t178 * MDP(17) + t170;
t168 = -0.2e1 * pkin(4);
t155 = qJ(5) * t162;
t154 = qJ(5) * t158;
t139 = t146 * t162;
t138 = t146 * t158;
t124 = t127 * t162;
t123 = t127 * t158;
t106 = pkin(4) * t119 + t171;
t100 = pkin(5) * t120 + t107;
t99 = t206 * t119 + t171;
t96 = t100 * t158 + t162 * t99;
t95 = t100 * t162 - t158 * t99;
t1 = [(t106 ^ 2 + t107 ^ 2 + t108 ^ 2) * MDP(28) + pkin(1) * MDP(9) * t212 + t135 * MDP(16) * t213 + MDP(1) + (MDP(24) * t125 - MDP(27) * t106) * t188 + (MDP(18) + MDP(33)) * t120 ^ 2 + (MDP(29) * t156 + 0.2e1 * t187) * t119 ^ 2 + (-0.2e1 * MDP(10) * pkin(1) + MDP(4) * t161 + MDP(5) * t212) * t161 + (MDP(11) * t136 - 0.2e1 * MDP(12) * t135 + MDP(17) * t213) * t136 + (0.2e1 * t125 * MDP(23) + t106 * t211 + (-MDP(19) + t177) * t188) * t119 + 0.2e1 * (t107 * t120 - t108 * t119) * MDP(25) + (-t119 * t98 + t120 * t95) * t208 + (t119 * t97 - t120 * t96) * t207; t161 * MDP(6) + (t120 * t129 - t199) * MDP(25) + (t107 * t129 + t108 * t127) * MDP(28) + (t180 * MDP(35) - t193) * t158 + t169 + (-t180 * t162 + t97) * MDP(34) + (-MDP(10) * t165 - MDP(9) * t161) * pkin(7) + t165 * MDP(7); MDP(8) + 0.2e1 * (-t160 * MDP(17) + t190) * pkin(2) - 0.2e1 * t192 - 0.2e1 * t191 + t183 + (t210 + t214) * t129 + (MDP(28) * t127 + t172) * t127; (t120 * t149 - t198) * MDP(25) + (t107 * t149 + t108 * t146) * MDP(28) + (t179 * MDP(35) - t193) * t158 + t169 + (-t179 * t162 + t97) * MDP(34); MDP(15) + t182 * MDP(23) + (t168 - t182) * MDP(26) + (t127 * t146 + t129 * t149) * MDP(28) + (t138 + t123) * MDP(34) + (t139 + t124) * MDP(35) + (t190 + (t189 * t163 - MDP(17)) * t160) * pkin(2) + t173 + (pkin(3) + t150) * t215; (t210 + t195) * t149 + 0.2e1 * (t163 * MDP(23) - t159 * MDP(24)) * pkin(3) + (MDP(28) * t146 + t172) * t146 + t183; (-pkin(4) * t120 - t200) * MDP(25) + (-pkin(4) * t107 + qJ(5) * t108) * MDP(28) + t170 + (t181 * MDP(35) - t193) * t158 + (-t181 * t162 + t97) * MDP(34); -t192 - t191 + (t168 + t196) * MDP(26) + (t167 + t174) * MDP(27) + (-pkin(4) * t129 + qJ(5) * t127) * MDP(28) + (t154 + t123) * MDP(34) + (t155 + t124) * MDP(35) + t186; t168 * MDP(26) + (-pkin(4) * t149 + qJ(5) * t146) * MDP(28) + (t154 + t138) * MDP(34) + (t155 + t139) * MDP(35) + (t216 * t163 + t215) * pkin(3) + t173; (t211 + t201) * pkin(4) + (MDP(28) * qJ(5) + t172) * qJ(5) + t186; t107 * MDP(28) + (MDP(25) + t176) * t120; MDP(26) + t214; MDP(26) + t195; MDP(26) - t201; MDP(28); t120 * MDP(33) + t95 * MDP(34) - t96 * MDP(35) + t177 * t119; t176 * t126 + t184; t176 * t145 + t184; -t206 * t194 + t153 + (MDP(35) * t206 - MDP(32)) * t158; t176; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
