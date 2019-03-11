% Calculate joint inertia matrix for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% MDP [34x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(34,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [34 1]), ...
  'S6RPRRRR9_inertiaJ_mdp_slag_vp: MDP has to be [34x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:25:31
% EndTime: 2019-03-09 07:25:33
% DurationCPUTime: 0.74s
% Computational Cost: add. (868->185), mult. (1669->253), div. (0->0), fcn. (1775->8), ass. (0->93)
t164 = sin(qJ(4));
t168 = cos(qJ(4));
t208 = pkin(8) + pkin(9);
t149 = t208 * t164;
t150 = t208 * t168;
t163 = sin(qJ(5));
t167 = cos(qJ(5));
t124 = -t167 * t149 - t163 * t150;
t125 = -t163 * t149 + t167 * t150;
t145 = t163 * t164 - t167 * t168;
t146 = t163 * t168 + t167 * t164;
t112 = -t146 * pkin(10) + t124;
t113 = -t145 * pkin(10) + t125;
t162 = sin(qJ(6));
t166 = cos(qJ(6));
t119 = t166 * t145 + t162 * t146;
t120 = -t162 * t145 + t166 * t146;
t179 = t120 * MDP(30) - t119 * MDP(31) + (t166 * t112 - t162 * t113) * MDP(33) - (t162 * t112 + t166 * t113) * MDP(34);
t172 = t146 * MDP(23) - t145 * MDP(24) + t124 * MDP(26) - t125 * MDP(27) + t179;
t214 = t164 * MDP(19) + t168 * MDP(20);
t217 = t164 * MDP(16) + t168 * MDP(17) - pkin(8) * t214 + t172;
t169 = cos(qJ(3));
t135 = t146 * t169;
t137 = t145 * t169;
t109 = t166 * t135 - t162 * t137;
t111 = -t162 * t135 - t166 * t137;
t216 = t111 * MDP(30) - t109 * MDP(31);
t215 = -t137 * MDP(23) - t135 * MDP(24);
t212 = -2 * MDP(22);
t211 = 0.2e1 * MDP(27);
t210 = -2 * MDP(29);
t209 = 0.2e1 * MDP(34);
t207 = (pkin(1) * MDP(6));
t206 = pkin(4) * t163;
t205 = pkin(9) * t169;
t165 = sin(qJ(3));
t204 = t165 * pkin(5);
t203 = t166 * pkin(5);
t202 = t167 * pkin(4);
t148 = t165 * pkin(3) - t169 * pkin(8) + qJ(2);
t143 = t168 * t148;
t170 = -pkin(1) - pkin(7);
t199 = t164 * t170;
t117 = -t168 * t205 + t143 + (pkin(4) - t199) * t165;
t197 = t168 * t170;
t183 = t165 * t197;
t121 = t183 + (t148 - t205) * t164;
t198 = t167 * t121;
t103 = t163 * t117 + t198;
t97 = -t135 * pkin(10) + t103;
t201 = t166 * t97;
t200 = t164 * t168;
t134 = t146 * t165;
t136 = t145 * t165;
t108 = -t166 * t134 + t162 * t136;
t110 = -t162 * t134 - t166 * t136;
t196 = t108 * MDP(33) - t110 * MDP(34);
t194 = t119 * MDP(33);
t193 = t120 * MDP(28);
t154 = pkin(5) + t202;
t151 = t166 * t154;
t138 = -t162 * t206 + t151;
t192 = t138 * MDP(33);
t139 = t162 * t154 + t166 * t206;
t191 = t139 * MDP(34);
t190 = t145 * MDP(26);
t189 = t146 * MDP(21);
t188 = t162 * MDP(34);
t186 = t167 * MDP(26);
t184 = MDP(25) + MDP(32);
t182 = t165 * MDP(32) + t216;
t155 = -t168 * pkin(4) - pkin(3);
t181 = MDP(15) * t200;
t102 = t167 * t117 - t163 * t121;
t96 = t137 * pkin(10) + t102 + t204;
t93 = -t162 * t97 + t166 * t96;
t180 = MDP(18) + t184;
t144 = (pkin(4) * t164 - t170) * t169;
t178 = -t134 * MDP(26) + t136 * MDP(27) + t196;
t94 = t162 * t96 + t201;
t177 = t168 * MDP(16) - t164 * MDP(17);
t176 = t168 * MDP(19) - t164 * MDP(20);
t174 = t165 * MDP(25) + t182 + t215;
t173 = (t166 * MDP(33) - t188) * pkin(5);
t161 = t169 ^ 2;
t160 = t168 ^ 2;
t159 = t165 ^ 2;
t158 = t164 ^ 2;
t129 = t145 * pkin(5) + t155;
t127 = t164 * t148 + t183;
t126 = -t165 * t199 + t143;
t118 = t135 * pkin(5) + t144;
t1 = [MDP(1) - (-t137 * MDP(21) + t135 * t212) * t137 + (t111 * MDP(28) + t109 * t210) * t111 + ((-2 * MDP(4) + t207) * pkin(1)) + (0.2e1 * t169 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + (t160 * MDP(14) + MDP(7) - 0.2e1 * t181) * t161 + t180 * t159 + 0.2e1 * (qJ(2) * MDP(12) + (-MDP(8) + t177) * t169 + t215 + t216) * t165 + 0.2e1 * (t118 * t109 + t93 * t165) * MDP(33) + (t118 * t111 - t94 * t165) * t209 + 0.2e1 * (t102 * t165 + t144 * t135) * MDP(26) + (-t103 * t165 - t144 * t137) * t211 + 0.2e1 * (t126 * t165 - t161 * t199) * MDP(19) + 0.2e1 * (-t127 * t165 - t161 * t197) * MDP(20); MDP(4) - t207 + (-t134 * t165 - t169 * t135) * MDP(26) + (t136 * t165 + t169 * t137) * MDP(27) + (t108 * t165 - t169 * t109) * MDP(33) + (-t110 * t165 - t169 * t111) * MDP(34) + t214 * (-t159 - t161); MDP(6); -t137 * t189 + (-t146 * t135 + t137 * t145) * MDP(22) + (t155 * t135 + t144 * t145) * MDP(26) + (-t155 * t137 + t144 * t146) * MDP(27) + t111 * t193 + (-t120 * t109 - t111 * t119) * MDP(29) + (t129 * t109 + t118 * t119) * MDP(33) + (t129 * t111 + t118 * t120) * MDP(34) + (MDP(9) + t170 * MDP(12) + MDP(14) * t200 + (-t158 + t160) * MDP(15) + (-pkin(3) * t164 + t197) * MDP(19) + (-pkin(3) * t168 - t199) * MDP(20)) * t169 + (-t170 * MDP(13) - MDP(10) + t217) * t165; -t165 * MDP(13) + (-t146 * MDP(27) - t120 * MDP(34) + MDP(12) + t176 - t190 - t194) * t169; 0.2e1 * t181 + 0.2e1 * t155 * t190 + 0.2e1 * t129 * t194 + t158 * MDP(14) + MDP(11) + 0.2e1 * t176 * pkin(3) + (t145 * t212 + t155 * t211 + t189) * t146 + (t119 * t210 + t129 * t209 + t193) * t120; t165 * MDP(18) + t126 * MDP(19) - t127 * MDP(20) + (t165 * t202 + t102) * MDP(26) + (-t198 + (-t165 * pkin(4) - t117) * t163) * MDP(27) + (t138 * t165 + t93) * MDP(33) + (-t139 * t165 - t94) * MDP(34) + t177 * t169 + t174; -t165 * t214 + t178; t217; 0.2e1 * (-t163 * MDP(27) + t186) * pkin(4) + 0.2e1 * t192 - 0.2e1 * t191 + t180; t102 * MDP(26) - t103 * MDP(27) + (t165 * t203 + t93) * MDP(33) + (-t201 + (-t96 - t204) * t162) * MDP(34) + t174; t178; t172; (t151 + t203) * MDP(33) + (-pkin(5) - t154) * t188 + (t186 + (-MDP(33) * t162 - MDP(34) * t166 - MDP(27)) * t163) * pkin(4) + t184; 0.2e1 * t173 + t184; t93 * MDP(33) - t94 * MDP(34) + t182; t196; t179; MDP(32) - t191 + t192; MDP(32) + t173; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
