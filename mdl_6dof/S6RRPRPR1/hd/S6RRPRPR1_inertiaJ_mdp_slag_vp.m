% Calculate joint inertia matrix for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR1_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:09:59
% EndTime: 2019-03-09 10:10:01
% DurationCPUTime: 0.66s
% Computational Cost: add. (1531->164), mult. (2830->230), div. (0->0), fcn. (3323->10), ass. (0->90)
t180 = cos(pkin(10));
t168 = t180 * pkin(2) + pkin(3);
t182 = sin(qJ(4));
t178 = sin(pkin(10));
t222 = t178 * pkin(2);
t223 = cos(qJ(4));
t149 = -t182 * t168 - t223 * t222;
t146 = qJ(5) - t149;
t177 = sin(pkin(11));
t179 = cos(pkin(11));
t211 = t177 ^ 2 + t179 ^ 2;
t214 = t211 * t146;
t181 = sin(qJ(6));
t184 = cos(qJ(6));
t156 = t181 * t177 - t184 * t179;
t158 = t184 * t177 + t181 * t179;
t213 = t158 * MDP(26) - t156 * MDP(27);
t191 = t156 * MDP(29) + t158 * MDP(30);
t228 = t179 * MDP(20) - t177 * MDP(21);
t190 = -t228 + t191;
t185 = cos(qJ(2));
t170 = -t185 * pkin(2) - pkin(1);
t183 = sin(qJ(2));
t196 = t178 * t183 - t180 * t185;
t144 = t196 * pkin(3) + t170;
t227 = 0.2e1 * t144;
t226 = 0.2e1 * t185;
t225 = 2 * MDP(22);
t224 = -2 * MDP(25);
t221 = t179 * pkin(5);
t174 = t179 * pkin(9);
t219 = -qJ(3) - pkin(7);
t218 = pkin(4) * MDP(23);
t164 = t219 * t183;
t165 = t219 * t185;
t139 = t180 * t164 + t178 * t165;
t157 = t178 * t185 + t180 * t183;
t128 = -t157 * pkin(8) + t139;
t140 = t178 * t164 - t180 * t165;
t129 = -t196 * pkin(8) + t140;
t120 = -t223 * t128 + t182 * t129;
t217 = t120 * t179;
t133 = t223 * t157 - t182 * t196;
t216 = t133 * t177;
t132 = t182 * t157 + t223 * t196;
t119 = t132 * pkin(4) - t133 * qJ(5) + t144;
t121 = t182 * t128 + t223 * t129;
t108 = t177 * t119 + t179 * t121;
t212 = t211 * qJ(5);
t123 = t156 * t133;
t210 = MDP(24) * t123;
t122 = t158 * t133;
t209 = t122 * MDP(27);
t208 = t123 * MDP(26);
t207 = t132 * MDP(28);
t148 = t223 * t168 - t182 * t222;
t147 = -pkin(4) - t148;
t206 = t147 * MDP(23);
t205 = t148 * MDP(18);
t204 = t149 * MDP(19);
t202 = MDP(17) + (MDP(24) * t158 + t156 * t224) * t158;
t107 = t179 * t119 - t177 * t121;
t201 = t211 * MDP(23);
t200 = -pkin(4) * t133 - qJ(5) * t132;
t199 = t107 * t179 + t108 * t177;
t198 = -t107 * t177 + t108 * t179;
t197 = -t132 * t146 + t133 * t147;
t195 = 0.2e1 * t228;
t194 = t177 * MDP(20) + t179 * MDP(21);
t104 = t132 * pkin(5) - t133 * t174 + t107;
t105 = -pkin(9) * t216 + t108;
t193 = (t184 * t104 - t181 * t105) * MDP(29) - (t181 * t104 + t184 * t105) * MDP(30);
t192 = t122 * MDP(29) - t123 * MDP(30);
t189 = 0.2e1 * t191;
t188 = t198 * MDP(22) + (-t158 * t122 + t123 * t156) * MDP(25) - t158 * t210 - t120 * MDP(18) - t121 * MDP(19) + t133 * MDP(15) + (-MDP(16) + t213) * t132;
t169 = -pkin(4) - t221;
t163 = t179 * qJ(5) + t174;
t162 = (-pkin(9) - qJ(5)) * t177;
t141 = t147 - t221;
t138 = t179 * t146 + t174;
t137 = (-pkin(9) - t146) * t177;
t136 = t181 * t162 + t184 * t163;
t135 = t184 * t162 - t181 * t163;
t125 = t181 * t137 + t184 * t138;
t124 = t184 * t137 - t181 * t138;
t116 = t120 * t177;
t111 = pkin(5) * t216 + t120;
t110 = t111 * t158;
t109 = t111 * t156;
t1 = [MDP(1) + pkin(1) * MDP(9) * t226 + (t139 ^ 2 + t140 ^ 2 + t170 ^ 2) * MDP(12) + (t107 ^ 2 + t108 ^ 2 + t120 ^ 2) * MDP(23) + (MDP(13) * t133 + MDP(19) * t227) * t133 - (t122 * t224 - t210) * t123 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t183 + MDP(5) * t226) * t183 + (-0.2e1 * t133 * MDP(14) + MDP(18) * t227 + t207 - 0.2e1 * t208 - 0.2e1 * t209) * t132 + 0.2e1 * (-t139 * t157 - t140 * t196) * MDP(11) + 0.2e1 * t192 * t111 + 0.2e1 * (t107 * MDP(20) - t108 * MDP(21) + t193) * t132 + 0.2e1 * (-t199 * MDP(22) + t194 * t120) * t133; t188 + (t197 * t177 - t217) * MDP(20) + (t120 * t147 + t198 * t146) * MDP(23) + (t197 * t179 + t116) * MDP(21) + t183 * MDP(6) + t185 * MDP(7) + (t141 * t122 + t124 * t132 + t109) * MDP(29) + (-t141 * t123 - t125 * t132 + t110) * MDP(30) + (-t185 * MDP(10) - t183 * MDP(9)) * pkin(7) + ((-t180 * t157 - t178 * t196) * MDP(11) + (t139 * t180 + t140 * t178) * MDP(12)) * pkin(2); MDP(8) + (t178 ^ 2 + t180 ^ 2) * MDP(12) * pkin(2) ^ 2 + t146 ^ 2 * t201 + (-t195 + t206) * t147 + t141 * t189 + 0.2e1 * t205 + 0.2e1 * t204 + t214 * t225 + t202; t170 * MDP(12) + t199 * MDP(23) + (-t211 * MDP(22) + MDP(19)) * t133 + (MDP(18) - t190) * t132; 0; MDP(12) + t201; (t200 * t177 - t217) * MDP(20) + (t200 * t179 + t116) * MDP(21) + (-t120 * pkin(4) + t198 * qJ(5)) * MDP(23) + (t169 * t122 + t135 * t132 + t109) * MDP(29) + (-t169 * t123 - t136 * t132 + t110) * MDP(30) + t188; t205 + t204 + (t212 + t214) * MDP(22) + (-t147 * pkin(4) + qJ(5) * t214) * MDP(23) + t202 + t228 * (pkin(4) - t147) + t191 * (t141 + t169); 0; t212 * t225 + qJ(5) ^ 2 * t201 + t169 * t189 + (t195 + t218) * pkin(4) + t202; t120 * MDP(23) + t194 * t133 + t192; t190 + t206; 0; t190 - t218; MDP(23); t193 + t207 - t208 - t209; t124 * MDP(29) - t125 * MDP(30) + t213; -t191; t135 * MDP(29) - t136 * MDP(30) + t213; 0; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
