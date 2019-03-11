% Calculate joint inertia matrix for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRR3_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:02:32
% EndTime: 2019-03-09 07:02:34
% DurationCPUTime: 0.67s
% Computational Cost: add. (898->160), mult. (1736->223), div. (0->0), fcn. (1837->10), ass. (0->87)
t161 = sin(qJ(4));
t165 = cos(qJ(4));
t205 = pkin(8) + pkin(9);
t144 = t205 * t161;
t145 = t205 * t165;
t160 = sin(qJ(5));
t164 = cos(qJ(5));
t122 = -t164 * t144 - t160 * t145;
t123 = -t160 * t144 + t164 * t145;
t140 = t160 * t161 - t164 * t165;
t141 = t160 * t165 + t164 * t161;
t110 = -t141 * pkin(10) + t122;
t111 = -t140 * pkin(10) + t123;
t159 = sin(qJ(6));
t163 = cos(qJ(6));
t116 = t163 * t140 + t159 * t141;
t117 = -t159 * t140 + t163 * t141;
t179 = t117 * MDP(28) - t116 * MDP(29) + (t163 * t110 - t159 * t111) * MDP(31) - (t159 * t110 + t163 * t111) * MDP(32);
t169 = t141 * MDP(21) - t140 * MDP(22) + t122 * MDP(24) - t123 * MDP(25) + t179;
t174 = t161 * MDP(17) + t165 * MDP(18);
t216 = t161 * MDP(14) + t165 * MDP(15) - t174 * pkin(8) + t169;
t186 = t159 * MDP(32);
t170 = (t163 * MDP(31) - t186) * pkin(5);
t184 = t164 * MDP(24);
t215 = (-t160 * MDP(25) + t184) * pkin(4);
t162 = sin(qJ(3));
t130 = t141 * t162;
t131 = t140 * t162;
t107 = t163 * t130 - t159 * t131;
t108 = -t159 * t130 - t163 * t131;
t195 = t108 * MDP(28) - t107 * MDP(29);
t214 = -t131 * MDP(21) - t130 * MDP(22);
t158 = cos(pkin(11));
t149 = -t158 * pkin(1) - pkin(2);
t166 = cos(qJ(3));
t139 = -t166 * pkin(3) - t162 * pkin(8) + t149;
t132 = t165 * t139;
t157 = sin(pkin(11));
t148 = t157 * pkin(1) + pkin(7);
t202 = t148 * t161;
t203 = pkin(9) * t162;
t109 = -t165 * t203 + t132 + (-pkin(4) - t202) * t166;
t200 = t148 * t166;
t182 = t165 * t200;
t112 = t182 + (t139 - t203) * t161;
t100 = t164 * t109 - t160 * t112;
t101 = t160 * t109 + t164 * t112;
t94 = -t166 * pkin(5) + t131 * pkin(10) + t100;
t95 = -t130 * pkin(10) + t101;
t91 = -t159 * t95 + t163 * t94;
t92 = t159 * t94 + t163 * t95;
t212 = t91 * MDP(31) - t92 * MDP(32) + t195;
t213 = t100 * MDP(24) - t101 * MDP(25) + t212 + t214;
t209 = -2 * MDP(20);
t208 = 0.2e1 * MDP(25);
t207 = -2 * MDP(27);
t206 = 0.2e1 * MDP(32);
t204 = pkin(4) * t160;
t201 = t148 * t165;
t199 = t161 * t165;
t196 = -t107 * MDP(31) - t108 * MDP(32);
t133 = (pkin(4) * t161 + t148) * t162;
t192 = t116 * MDP(31);
t191 = t117 * MDP(26);
t151 = t164 * pkin(4) + pkin(5);
t146 = t163 * t151;
t190 = (-t159 * t204 + t146) * MDP(31);
t189 = (t159 * t151 + t163 * t204) * MDP(32);
t188 = t140 * MDP(24);
t187 = t141 * MDP(19);
t185 = t162 * MDP(11);
t183 = MDP(23) + MDP(30);
t152 = -t165 * pkin(4) - pkin(3);
t181 = MDP(13) * t199;
t180 = MDP(16) + t183;
t178 = -t130 * MDP(24) + t131 * MDP(25) + t196;
t177 = t165 * MDP(14) - t161 * MDP(15);
t175 = t165 * MDP(17) - t161 * MDP(18);
t171 = MDP(30) - t189 + t190;
t155 = t165 ^ 2;
t154 = t162 ^ 2;
t153 = t161 ^ 2;
t125 = t140 * pkin(5) + t152;
t119 = t161 * t139 + t182;
t118 = -t161 * t200 + t132;
t113 = t130 * pkin(5) + t133;
t1 = [0.2e1 * t149 * t185 + MDP(1) + (t157 ^ 2 + t158 ^ 2) * MDP(4) * pkin(1) ^ 2 - (-t131 * MDP(19) + t130 * t209) * t131 + (t108 * MDP(26) + t107 * t207) * t108 + t180 * t166 ^ 2 + (t155 * MDP(12) + MDP(5) - 0.2e1 * t181) * t154 + 0.2e1 * (-t149 * MDP(10) + (MDP(6) - t177) * t162 - t214 - t195) * t166 + 0.2e1 * (t113 * t107 - t91 * t166) * MDP(31) + (t113 * t108 + t92 * t166) * t206 + 0.2e1 * (-t100 * t166 + t133 * t130) * MDP(24) + (t101 * t166 - t133 * t131) * t208 + 0.2e1 * (-t118 * t166 + t154 * t202) * MDP(17) + 0.2e1 * (t119 * t166 + t154 * t201) * MDP(18); 0; MDP(4); -t131 * t187 + (-t141 * t130 + t131 * t140) * MDP(20) + (t152 * t130 + t133 * t140) * MDP(24) + (-t152 * t131 + t133 * t141) * MDP(25) + t108 * t191 + (-t117 * t107 - t108 * t116) * MDP(27) + (t125 * t107 + t113 * t116) * MDP(31) + (t125 * t108 + t113 * t117) * MDP(32) + (MDP(7) - t148 * MDP(10) + MDP(12) * t199 + (-t153 + t155) * MDP(13) + (-pkin(3) * t161 - t201) * MDP(17) + (-pkin(3) * t165 + t202) * MDP(18)) * t162 + (-t148 * MDP(11) + MDP(8) - t216) * t166; -t185 + (-t141 * MDP(25) - t117 * MDP(32) + MDP(10) + t175 - t188 - t192) * t166; 0.2e1 * t181 + 0.2e1 * t152 * t188 + 0.2e1 * t125 * t192 + t153 * MDP(12) + MDP(9) + 0.2e1 * t175 * pkin(3) + (t140 * t209 + t152 * t208 + t187) * t141 + (t116 * t207 + t125 * t206 + t191) * t117; t118 * MDP(17) - t119 * MDP(18) + t177 * t162 + (-MDP(16) - MDP(23) - t171 - t215) * t166 + t213; -t174 * t162 + t178; t216; t180 - 0.2e1 * t189 + 0.2e1 * t190 + 0.2e1 * t215; (-t183 - t170) * t166 + t213; t178; t169; (t163 * pkin(5) + t146) * MDP(31) + (-pkin(5) - t151) * t186 + (t184 + (-MDP(31) * t159 - MDP(32) * t163 - MDP(25)) * t160) * pkin(4) + t183; 0.2e1 * t170 + t183; -t166 * MDP(30) + t212; t196; t179; t171; MDP(30) + t170; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
