% Calculate joint inertia matrix for
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR1_inertiaJ_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:29:51
% EndTime: 2019-03-10 03:29:53
% DurationCPUTime: 0.74s
% Computational Cost: add. (1814->194), mult. (3314->240), div. (0->0), fcn. (3934->10), ass. (0->100)
t171 = sin(qJ(6));
t176 = cos(qJ(6));
t215 = t171 * MDP(34) + t176 * MDP(35);
t214 = MDP(37) * t176;
t183 = -0.2e1 * MDP(38) * t171 + 0.2e1 * t214;
t174 = sin(qJ(3));
t175 = sin(qJ(2));
t179 = cos(qJ(2));
t222 = cos(qJ(3));
t148 = t174 * t175 - t222 * t179;
t203 = -t179 * pkin(2) - pkin(1);
t138 = t148 * pkin(3) + t203;
t149 = t174 * t179 + t222 * t175;
t173 = sin(qJ(4));
t178 = cos(qJ(4));
t190 = t178 * t148 + t173 * t149;
t121 = t190 * pkin(4) + t138;
t226 = 0.2e1 * t121;
t225 = 0.2e1 * t138;
t224 = 0.2e1 * t179;
t223 = pkin(7) + pkin(8);
t221 = pkin(2) * t174;
t220 = pkin(5) * t171;
t219 = t173 * pkin(3);
t133 = -t173 * t148 + t178 * t149;
t151 = t223 * t175;
t152 = t223 * t179;
t197 = -t222 * t151 - t174 * t152;
t128 = -t149 * pkin(9) + t197;
t185 = t174 * t151 - t222 * t152;
t129 = -t148 * pkin(9) - t185;
t199 = t178 * t128 - t173 * t129;
t109 = -t133 * pkin(10) + t199;
t191 = -t173 * t128 - t178 * t129;
t110 = -t190 * pkin(10) - t191;
t172 = sin(qJ(5));
t177 = cos(qJ(5));
t105 = -t177 * t109 + t172 * t110;
t102 = t105 * t171;
t218 = t105 * t176;
t217 = t171 * t176;
t161 = t222 * pkin(2) + pkin(3);
t145 = t173 * t161 + t178 * t221;
t216 = t177 * t145;
t119 = t172 * t133 + t177 * t190;
t213 = t119 * MDP(36);
t154 = t178 * t161;
t143 = -t173 * t221 + t154;
t140 = pkin(4) + t143;
t136 = t177 * t140;
t198 = t172 * t145 - t136;
t212 = t198 * MDP(30);
t127 = -t172 * t140 - t216;
t211 = t127 * MDP(31);
t167 = t178 * pkin(3);
t160 = t167 + pkin(4);
t153 = t177 * t160;
t142 = -t172 * t219 + t153;
t210 = t142 * MDP(30);
t209 = t143 * MDP(23);
t144 = -t172 * t160 - t177 * t219;
t208 = t144 * MDP(31);
t207 = t145 * MDP(24);
t206 = t172 * MDP(31);
t205 = t177 * MDP(31);
t204 = t178 * MDP(23);
t202 = MDP(33) * t217;
t169 = t171 ^ 2;
t201 = t169 * MDP(32) + MDP(29) + 0.2e1 * t202;
t200 = t222 * MDP(16);
t196 = MDP(22) + t201;
t120 = t177 * t133 - t172 * t190;
t195 = -pkin(5) * t120 - pkin(11) * t119;
t124 = -pkin(5) + t198;
t125 = pkin(11) - t127;
t194 = -t119 * t125 + t120 * t124;
t139 = -pkin(5) - t142;
t141 = pkin(11) - t144;
t193 = -t119 * t141 + t120 * t139;
t158 = t172 * pkin(4) + pkin(11);
t166 = t177 * pkin(4);
t159 = -t166 - pkin(5);
t192 = -t119 * t158 + t120 * t159;
t189 = MDP(15) + t196;
t188 = MDP(34) * t176 - MDP(35) * t171;
t186 = -MDP(37) * t171 - MDP(38) * t176;
t184 = (t177 * MDP(30) - t206) * pkin(4);
t106 = t172 * t109 + t177 * t110;
t170 = t176 ^ 2;
t182 = -t105 * MDP(30) - t106 * MDP(31) + ((-t169 + t170) * MDP(33) + MDP(32) * t217 + MDP(27)) * t120 + (-MDP(28) + t215) * t119;
t181 = t133 * MDP(20) - t190 * MDP(21) + t199 * MDP(23) + t191 * MDP(24) + t182;
t180 = t149 * MDP(13) - t148 * MDP(14) + t197 * MDP(16) + t185 * MDP(17) + t181;
t168 = pkin(5) * t176;
t155 = t159 * t171;
t137 = t139 * t171;
t122 = t124 * t171;
t107 = t119 * pkin(5) - t120 * pkin(11) + t121;
t101 = t176 * t106 + t171 * t107;
t100 = -t171 * t106 + t176 * t107;
t1 = [pkin(1) * MDP(9) * t224 + t190 * MDP(23) * t225 + t120 * MDP(31) * t226 + MDP(1) + (MDP(30) * t226 + t213) * t119 + 0.2e1 * (t100 * t119 + t120 * t102) * MDP(37) + 0.2e1 * (-t101 * t119 + t120 * t218) * MDP(38) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t175 + MDP(5) * t224) * t175 + (MDP(18) * t133 - 0.2e1 * t190 * MDP(19) + MDP(24) * t225) * t133 + (t170 * MDP(32) + MDP(25) - 0.2e1 * t202) * t120 ^ 2 + (MDP(11) * t149 - 0.2e1 * t148 * MDP(12)) * t149 + 0.2e1 * (t148 * MDP(16) + t149 * MDP(17)) * t203 + 0.2e1 * (-MDP(26) + t188) * t120 * t119; (t194 * t171 - t218) * MDP(37) + (t194 * t176 + t102) * MDP(38) + t180 + t179 * MDP(7) + t175 * MDP(6) + (-t179 * MDP(10) - t175 * MDP(9)) * pkin(7); MDP(8) - t124 * t183 + 0.2e1 * (-t174 * MDP(17) + t200) * pkin(2) + 0.2e1 * t209 - 0.2e1 * t207 - 0.2e1 * t212 + 0.2e1 * t211 + t189; (t193 * t171 - t218) * MDP(37) + (t193 * t176 + t102) * MDP(38) + t180; (t154 + t167) * MDP(23) + (t136 + t153) * MDP(30) - t145 * t205 + (t137 + t122) * MDP(38) + (-t124 - t139) * t214 + ((-pkin(3) - t161) * MDP(24) - pkin(3) * t205) * t173 + ((-t145 - t219) * MDP(30) + (-t140 - t160) * MDP(31)) * t172 + (t200 + (-MDP(23) * t173 - MDP(24) * t178 - MDP(17)) * t174) * pkin(2) + t189; -t139 * t183 + 0.2e1 * (-t173 * MDP(24) + t204) * pkin(3) + 0.2e1 * t210 + 0.2e1 * t208 + t189; (t192 * t171 - t218) * MDP(37) + (t192 * t176 + t102) * MDP(38) + t181; t209 - t207 + (t166 - t198) * MDP(30) + (-t216 + (-pkin(4) - t140) * t172) * MDP(31) + (t155 + t122) * MDP(38) + (-t124 - t159) * t214 + t196; (t153 + t166) * MDP(30) + (t155 + t137) * MDP(38) + (-t139 - t159) * t214 + (-pkin(4) - t160) * t206 + (t204 + (-t172 * MDP(30) - MDP(24) - t205) * t173) * pkin(3) + t196; -t159 * t183 + 0.2e1 * t184 + t196; (t195 * t171 - t218) * MDP(37) + (t195 * t176 + t102) * MDP(38) + t182; -t212 + t211 + (-t124 * t176 + t168) * MDP(37) + (t122 - t220) * MDP(38) + t201; t210 + t208 + (-t139 * t176 + t168) * MDP(37) + (t137 - t220) * MDP(38) + t201; (-t159 * t176 + t168) * MDP(37) + (t155 - t220) * MDP(38) + t184 + t201; pkin(5) * t183 + t201; t100 * MDP(37) - t101 * MDP(38) + t188 * t120 + t213; t186 * t125 + t215; t186 * t141 + t215; t186 * t158 + t215; t186 * pkin(11) + t215; MDP(36);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
