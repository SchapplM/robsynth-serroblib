% Calculate joint inertia matrix for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRPR4_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:20:40
% EndTime: 2019-03-08 23:20:41
% DurationCPUTime: 0.66s
% Computational Cost: add. (869->173), mult. (1871->268), div. (0->0), fcn. (2085->12), ass. (0->85)
t163 = sin(qJ(4));
t167 = cos(qJ(4));
t195 = -qJ(5) - pkin(9);
t147 = t195 * t163;
t148 = t195 * t167;
t158 = sin(pkin(12));
t160 = cos(pkin(12));
t126 = t160 * t147 + t148 * t158;
t140 = t158 * t167 + t160 * t163;
t112 = -pkin(10) * t140 + t126;
t127 = t158 * t147 - t160 * t148;
t139 = -t158 * t163 + t160 * t167;
t113 = pkin(10) * t139 + t127;
t162 = sin(qJ(6));
t166 = cos(qJ(6));
t118 = -t166 * t139 + t140 * t162;
t119 = t139 * t162 + t140 * t166;
t177 = t119 * MDP(23) - t118 * MDP(24) + (t112 * t166 - t113 * t162) * MDP(26) - (t112 * t162 + t113 * t166) * MDP(27);
t205 = t177 - (MDP(17) * t163 + MDP(18) * t167) * pkin(9) + t163 * MDP(14) + t167 * MDP(15);
t164 = sin(qJ(3));
t132 = t140 * t164;
t189 = t164 * t167;
t191 = t163 * t164;
t133 = -t158 * t191 + t160 * t189;
t110 = t166 * t132 + t133 * t162;
t111 = -t132 * t162 + t133 * t166;
t186 = t111 * MDP(23) - t110 * MDP(24);
t168 = cos(qJ(3));
t146 = -pkin(3) * t168 - pkin(9) * t164 - pkin(2);
t141 = t167 * t146;
t198 = pkin(8) * t163;
t120 = -qJ(5) * t189 + t141 + (-pkin(4) - t198) * t168;
t196 = pkin(8) * t168;
t179 = t167 * t196;
t125 = t179 + (-qJ(5) * t164 + t146) * t163;
t104 = t160 * t120 - t125 * t158;
t100 = -pkin(5) * t168 - pkin(10) * t133 + t104;
t105 = t158 * t120 + t160 * t125;
t103 = -pkin(10) * t132 + t105;
t91 = t166 * t100 - t103 * t162;
t92 = t100 * t162 + t103 * t166;
t204 = t91 * MDP(26) - t92 * MDP(27) + t186;
t202 = 2 * MDP(19);
t201 = -2 * MDP(22);
t200 = 0.2e1 * MDP(27);
t199 = pkin(4) * t158;
t197 = pkin(8) * t167;
t161 = cos(pkin(6));
t159 = sin(pkin(6));
t165 = sin(qJ(2));
t193 = t159 * t165;
t138 = t161 * t164 + t168 * t193;
t169 = cos(qJ(2));
t192 = t159 * t169;
t123 = -t138 * t163 - t167 * t192;
t124 = t138 * t167 - t163 * t192;
t106 = t123 * t160 - t124 * t158;
t107 = t123 * t158 + t124 * t160;
t95 = t106 * t166 - t107 * t162;
t96 = t106 * t162 + t107 * t166;
t194 = t95 * MDP(26) - t96 * MDP(27);
t190 = t163 * t167;
t145 = pkin(4) * t191 + t164 * pkin(8);
t185 = MDP(11) * t164;
t184 = MDP(21) * t119;
t183 = MDP(26) * t118;
t151 = pkin(4) * t160 + pkin(5);
t182 = MDP(26) * (t151 * t166 - t162 * t199);
t181 = (t151 * t162 + t166 * t199) * MDP(27);
t180 = MDP(16) + MDP(25);
t152 = -pkin(4) * t167 - pkin(3);
t178 = MDP(13) * t190;
t176 = MDP(14) * t167 - MDP(15) * t163;
t174 = MDP(17) * t167 - MDP(18) * t163;
t172 = MDP(25) - t181 + t182;
t171 = MDP(20) * t152 + MDP(27) * t119 + t183;
t156 = t167 ^ 2;
t155 = t164 ^ 2;
t154 = t163 ^ 2;
t137 = -t161 * t168 + t164 * t193;
t131 = -pkin(5) * t139 + t152;
t130 = t146 * t163 + t179;
t129 = -t163 * t196 + t141;
t121 = pkin(5) * t132 + t145;
t1 = [MDP(1) + (t106 ^ 2 + t107 ^ 2 + t137 ^ 2) * MDP(20); (-t123 * t168 + t137 * t191) * MDP(17) + (t124 * t168 + t137 * t189) * MDP(18) + (-t106 * t133 - t107 * t132) * MDP(19) + (t104 * t106 + t105 * t107 + t137 * t145) * MDP(20) + (t110 * t137 - t168 * t95) * MDP(26) + (t111 * t137 + t168 * t96) * MDP(27) + (-t165 * MDP(4) + (MDP(10) * t168 + MDP(3) - t185) * t169) * t159; MDP(2) - 0.2e1 * pkin(2) * t185 + (t104 ^ 2 + t105 ^ 2 + t145 ^ 2) * MDP(20) + t180 * t168 ^ 2 + (MDP(21) * t111 + t110 * t201) * t111 + (MDP(12) * t156 + MDP(5) - 0.2e1 * t178) * t155 + 0.2e1 * (pkin(2) * MDP(10) + (MDP(6) - t176) * t164 - t186) * t168 + 0.2e1 * (-t129 * t168 + t155 * t198) * MDP(17) + 0.2e1 * (t130 * t168 + t155 * t197) * MDP(18) + (-t104 * t133 - t105 * t132) * t202 + 0.2e1 * (t121 * t110 - t91 * t168) * MDP(26) + (t111 * t121 + t92 * t168) * t200; -t138 * MDP(11) + (-t106 * t140 + t107 * t139) * MDP(19) + (t106 * t126 + t107 * t127) * MDP(20) + (-MDP(10) + t171 - t174) * t137; (-t104 * t140 + t105 * t139 - t126 * t133 - t127 * t132) * MDP(19) + (t104 * t126 + t105 * t127 + t145 * t152) * MDP(20) + t111 * t184 + (-t110 * t119 - t111 * t118) * MDP(22) + (t110 * t131 + t118 * t121) * MDP(26) + (t111 * t131 + t119 * t121) * MDP(27) + (-pkin(8) * MDP(11) + MDP(8) - t205) * t168 + (MDP(7) - pkin(8) * MDP(10) + MDP(12) * t190 + (-t154 + t156) * MDP(13) + (-pkin(3) * t163 - t197) * MDP(17) + (-pkin(3) * t167 + t198) * MDP(18)) * t164; MDP(9) + t154 * MDP(12) + 0.2e1 * t178 + (-t126 * t140 + t127 * t139) * t202 + (t126 ^ 2 + t127 ^ 2 + t152 ^ 2) * MDP(20) + 0.2e1 * t131 * t183 + 0.2e1 * t174 * pkin(3) + (t118 * t201 + t131 * t200 + t184) * t119; t123 * MDP(17) - t124 * MDP(18) + (t106 * t160 + t107 * t158) * MDP(20) * pkin(4) + t194; t129 * MDP(17) - t130 * MDP(18) + (-MDP(16) - t172) * t168 + t176 * t164 + ((-t132 * t158 - t133 * t160) * MDP(19) + (t104 * t160 + t105 * t158) * MDP(20)) * pkin(4) + t204; ((t139 * t158 - t140 * t160) * MDP(19) + (t126 * t160 + t127 * t158) * MDP(20)) * pkin(4) + t205; (t158 ^ 2 + t160 ^ 2) * MDP(20) * pkin(4) ^ 2 + 0.2e1 * t182 - 0.2e1 * t181 + t180; t137 * MDP(20); MDP(20) * t145 + t110 * MDP(26) + t111 * MDP(27); t171; 0; MDP(20); t194; -t168 * MDP(25) + t204; t177; t172; 0; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
