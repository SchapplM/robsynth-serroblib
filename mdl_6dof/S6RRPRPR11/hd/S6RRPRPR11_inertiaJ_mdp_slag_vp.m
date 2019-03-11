% Calculate joint inertia matrix for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR11_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:15:05
% EndTime: 2019-03-09 11:15:07
% DurationCPUTime: 0.61s
% Computational Cost: add. (914->167), mult. (1616->238), div. (0->0), fcn. (1668->8), ass. (0->84)
t226 = pkin(3) + pkin(7);
t180 = sin(qJ(4));
t183 = cos(qJ(4));
t185 = -pkin(2) - pkin(8);
t215 = -qJ(5) + t185;
t158 = t215 * t180;
t159 = t215 * t183;
t177 = sin(pkin(10));
t178 = cos(pkin(10));
t139 = -t158 * t177 + t178 * t159;
t192 = t177 * t180 - t178 * t183;
t121 = pkin(9) * t192 + t139;
t140 = t178 * t158 + t177 * t159;
t193 = t177 * t183 + t178 * t180;
t122 = -pkin(9) * t193 + t140;
t179 = sin(qJ(6));
t182 = cos(qJ(6));
t194 = -t179 * t192 + t182 * t193;
t199 = -t179 * t193 - t182 * t192;
t198 = t199 * MDP(26) - t194 * MDP(27) + (t121 * t182 - t122 * t179) * MDP(29) - (t121 * t179 + t122 * t182) * MDP(30);
t225 = (MDP(20) * t185 + MDP(17)) * t183 + (-MDP(21) * t185 - MDP(18)) * t180 + t198;
t184 = cos(qJ(2));
t143 = t192 * t184;
t144 = t193 * t184;
t123 = -t182 * t143 - t144 * t179;
t124 = t143 * t179 - t144 * t182;
t224 = t124 * MDP(26) - t123 * MDP(27);
t223 = 0.2e1 * t184;
t222 = -2 * MDP(25);
t221 = 0.2e1 * MDP(30);
t220 = pkin(4) * t177;
t219 = pkin(2) * MDP(14);
t181 = sin(qJ(2));
t161 = t226 * t181;
t218 = t161 * t180;
t217 = t180 * t183;
t216 = t183 * t184;
t168 = t180 * pkin(4) + qJ(3);
t157 = t183 * t161;
t202 = -qJ(3) * t181 - pkin(1);
t153 = t184 * t185 + t202;
t200 = qJ(5) * t184 - t153;
t125 = pkin(4) * t181 + t180 * t200 + t157;
t127 = -t183 * t200 + t218;
t116 = t177 * t125 + t178 * t127;
t214 = MDP(29) * t199 - MDP(30) * t194;
t162 = t226 * t184;
t167 = pkin(4) * t178 + pkin(5);
t145 = t167 * t182 - t179 * t220;
t213 = MDP(29) * t145;
t212 = t194 * MDP(29);
t211 = t199 * MDP(24);
t146 = t167 * t179 + t182 * t220;
t210 = t146 * MDP(30);
t209 = t180 * MDP(20);
t208 = t183 * MDP(21);
t207 = pkin(7) ^ 2 * MDP(14);
t206 = MDP(19) + MDP(28);
t205 = t181 * MDP(28) + t224;
t149 = pkin(4) * t216 + t162;
t204 = MDP(16) * t217;
t203 = t192 ^ 2 + t193 ^ 2;
t201 = MDP(12) - t219;
t115 = t178 * t125 - t127 * t177;
t109 = pkin(5) * t181 + pkin(9) * t144 + t115;
t110 = pkin(9) * t143 + t116;
t106 = t182 * t109 - t110 * t179;
t107 = t109 * t179 + t110 * t182;
t197 = -t115 * t192 + t116 * t193;
t196 = -t139 * t192 + t140 * t193;
t195 = t177 * t193 - t178 * t192;
t191 = -MDP(17) * t180 - MDP(18) * t183;
t190 = t183 * MDP(20) - t180 * MDP(21);
t189 = pkin(7) * MDP(14) + MDP(11) + t190;
t176 = t184 ^ 2;
t175 = t183 ^ 2;
t174 = t181 ^ 2;
t173 = t180 ^ 2;
t160 = -pkin(2) * t184 + t202;
t142 = pkin(5) * t193 + t168;
t138 = t153 * t183 + t218;
t137 = -t153 * t180 + t157;
t128 = -pkin(5) * t143 + t149;
t1 = [(t115 ^ 2 + t116 ^ 2 + t149 ^ 2) * MDP(23) + pkin(1) * MDP(9) * t223 + MDP(1) + (MDP(12) * t223 + MDP(14) * t160) * t160 + (t173 * MDP(15) + 0.2e1 * t204 + t207) * t176 + (t124 * MDP(24) + t123 * t222) * t124 + (MDP(4) + t206 + t207) * t174 + 0.2e1 * (-pkin(1) * MDP(10) - t160 * MDP(13) + (MDP(5) + t191) * t184 + t224) * t181 + 0.2e1 * (t115 * t144 + t116 * t143) * MDP(22) + 0.2e1 * (t106 * t181 + t123 * t128) * MDP(29) + (-t107 * t181 + t124 * t128) * t221 + 0.2e1 * (-t162 * t180 * t184 - t138 * t181) * MDP(21) + 0.2e1 * (t137 * t181 + t162 * t216) * MDP(20) + 0.2e1 * (t174 + t176) * MDP(11) * pkin(7); (t139 * t144 + t140 * t143 - t197) * MDP(22) + (t115 * t139 + t116 * t140 + t149 * t168) * MDP(23) + t124 * t211 + (-t123 * t199 - t124 * t194) * MDP(25) + (t123 * t142 + t128 * t194) * MDP(29) + (t124 * t142 + t128 * t199) * MDP(30) + (t208 + t209) * t162 + (MDP(7) - MDP(15) * t217 + (t173 - t175) * MDP(16) + (-MDP(10) + MDP(13)) * pkin(7) + t189 * qJ(3)) * t184 + (-pkin(2) * MDP(11) + MDP(6) + (-MDP(9) + t201) * pkin(7) + t225) * t181; MDP(8) + t175 * MDP(15) - 0.2e1 * t204 - 0.2e1 * t196 * MDP(22) + (t139 ^ 2 + t140 ^ 2 + t168 ^ 2) * MDP(23) + 0.2e1 * t142 * t212 + (-0.2e1 * MDP(12) + t219) * pkin(2) + (t142 * t221 + t194 * t222 + t211) * t199 + (qJ(3) * MDP(14) + 0.2e1 * MDP(13) + 0.2e1 * t208 + 0.2e1 * t209) * qJ(3); (t143 * t193 - t144 * t192) * MDP(22) + t197 * MDP(23) + (t189 + t214) * t181; -MDP(22) * t203 + MDP(23) * t196 + t201; MDP(23) * t203 + MDP(14); t181 * MDP(19) + t137 * MDP(20) - t138 * MDP(21) + (t145 * t181 + t106) * MDP(29) + (-t146 * t181 - t107) * MDP(30) + t191 * t184 + ((t143 * t177 + t144 * t178) * MDP(22) + (t115 * t178 + t116 * t177) * MDP(23)) * pkin(4) + t205; (-t195 * MDP(22) + (t139 * t178 + t140 * t177) * MDP(23)) * pkin(4) + t225; MDP(23) * pkin(4) * t195 + t190 + t214; (t177 ^ 2 + t178 ^ 2) * MDP(23) * pkin(4) ^ 2 + 0.2e1 * t213 - 0.2e1 * t210 + t206; MDP(23) * t149 + t123 * MDP(29) + t124 * MDP(30); MDP(23) * t168 + MDP(30) * t199 + t212; 0; 0; MDP(23); t106 * MDP(29) - t107 * MDP(30) + t205; t198; t214; MDP(28) - t210 + t213; 0; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
