% Calculate joint inertia matrix for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRPRR2_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:00:18
% EndTime: 2019-03-08 22:00:19
% DurationCPUTime: 0.48s
% Computational Cost: add. (676->142), mult. (1396->227), div. (0->0), fcn. (1618->12), ass. (0->77)
t181 = qJ(4) + pkin(8);
t137 = sin(pkin(12));
t139 = cos(pkin(12));
t142 = sin(qJ(3));
t176 = cos(qJ(3));
t124 = t137 * t176 + t139 * t142;
t140 = sin(qJ(6));
t141 = sin(qJ(5));
t144 = cos(qJ(6));
t145 = cos(qJ(5));
t127 = t140 * t145 + t144 * t141;
t101 = t127 * t124;
t126 = t140 * t141 - t144 * t145;
t102 = t126 * t124;
t180 = -t102 * MDP(23) - t101 * MDP(24);
t131 = t137 * pkin(3) + pkin(9);
t153 = t141 * MDP(19) + t145 * MDP(20);
t174 = pkin(10) + t131;
t119 = t174 * t141;
t120 = t174 * t145;
t155 = t127 * MDP(23) - t126 * MDP(24) + (-t144 * t119 - t140 * t120) * MDP(26) - (-t140 * t119 + t144 * t120) * MDP(27);
t179 = t141 * MDP(16) + t145 * MDP(17) - t153 * t131 + t155;
t154 = t145 * MDP(19) - t141 * MDP(20);
t115 = t126 * MDP(26);
t163 = -t127 * MDP(27) - t115;
t149 = t154 + t163;
t178 = -2 * MDP(22);
t177 = 0.2e1 * MDP(27);
t123 = t137 * t142 - t139 * t176;
t175 = t123 * pkin(5);
t138 = sin(pkin(6));
t143 = sin(qJ(2));
t167 = t138 * t143;
t170 = cos(pkin(6));
t113 = t170 * t142 + t176 * t167;
t150 = -t142 * t167 + t170 * t176;
t100 = t139 * t113 + t137 * t150;
t146 = cos(qJ(2));
t166 = t138 * t146;
t93 = -t141 * t100 - t145 * t166;
t94 = t145 * t100 - t141 * t166;
t86 = -t140 * t94 + t144 * t93;
t87 = t140 * t93 + t144 * t94;
t173 = t86 * MDP(26) - t87 * MDP(27);
t172 = MDP(13) * pkin(3);
t134 = -t176 * pkin(3) - pkin(2);
t108 = t123 * pkin(4) - t124 * pkin(9) + t134;
t129 = t181 * t176;
t156 = t181 * t142;
t111 = t139 * t129 - t137 * t156;
t164 = t145 * t111;
t90 = t164 + (-pkin(10) * t124 + t108) * t141;
t171 = t144 * t90;
t169 = t124 * t141;
t168 = t124 * t145;
t165 = t141 * t145;
t162 = t102 * MDP(21);
t161 = t134 * MDP(13);
t160 = MDP(18) + MDP(25);
t159 = t123 * MDP(25) + t180;
t132 = -t139 * pkin(3) - pkin(4);
t158 = MDP(15) * t165;
t157 = MDP(10) * t176;
t91 = t145 * t108 - t141 * t111;
t89 = -pkin(10) * t168 + t175 + t91;
t82 = -t140 * t90 + t144 * t89;
t109 = t137 * t129 + t139 * t156;
t152 = (MDP(26) * t144 - MDP(27) * t140) * pkin(5);
t151 = (t145 * MDP(16) - t141 * MDP(17)) * t124;
t136 = t145 ^ 2;
t135 = t141 ^ 2;
t128 = -t145 * pkin(5) + t132;
t98 = t137 * t113 - t139 * t150;
t95 = pkin(5) * t169 + t109;
t92 = t141 * t108 + t164;
t83 = t140 * t89 + t171;
t1 = [MDP(1) + (t138 ^ 2 * t146 ^ 2 + t100 ^ 2 + t98 ^ 2) * MDP(13); (-t100 * t123 + t98 * t124) * MDP(12) + (t100 * t111 + t98 * t109) * MDP(13) + (t93 * t123 + t98 * t169) * MDP(19) + (-t94 * t123 + t98 * t168) * MDP(20) + (t98 * t101 + t86 * t123) * MDP(26) + (-t98 * t102 - t87 * t123) * MDP(27) + (-t143 * MDP(4) + (-MDP(11) * t142 + MDP(3) + t157 - t161) * t146) * t138; MDP(2) + 0.2e1 * pkin(2) * t157 + (t109 ^ 2 + t111 ^ 2 + t134 ^ 2) * MDP(13) + (t136 * MDP(14) - 0.2e1 * t158) * t124 ^ 2 + t160 * t123 ^ 2 - (t101 * t178 - t162) * t102 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t142 + 0.2e1 * t176 * MDP(6)) * t142 + 0.2e1 * (t151 + t180) * t123 + 0.2e1 * (t109 * t124 - t111 * t123) * MDP(12) + 0.2e1 * (t109 * t169 + t91 * t123) * MDP(19) + 0.2e1 * (t109 * t168 - t92 * t123) * MDP(20) + 0.2e1 * (t95 * t101 + t82 * t123) * MDP(26) + (-t95 * t102 - t83 * t123) * t177; t150 * MDP(10) - t113 * MDP(11) + t100 * t137 * t172 + (-t139 * t172 - t149) * t98; t142 * MDP(7) + t176 * MDP(8) - t127 * t162 + (-t127 * t101 + t102 * t126) * MDP(22) + (t128 * t101 + t95 * t126) * MDP(26) + (-t128 * t102 + t95 * t127) * MDP(27) - t154 * t109 + (-t142 * MDP(10) - t176 * MDP(11)) * pkin(8) + (MDP(14) * t165 + (-t135 + t136) * MDP(15) + t153 * t132) * t124 + t179 * t123 + ((-t123 * t137 - t124 * t139) * MDP(12) + (-t109 * t139 + t111 * t137) * MDP(13)) * pkin(3); 0.2e1 * t158 + 0.2e1 * t128 * t115 + t135 * MDP(14) + MDP(9) + (t137 ^ 2 + t139 ^ 2) * MDP(13) * pkin(3) ^ 2 - 0.2e1 * t154 * t132 + (MDP(21) * t127 + t126 * t178 + t128 * t177) * t127; -MDP(13) * t166; t149 * t123 + t161; 0; MDP(13); t93 * MDP(19) - t94 * MDP(20) + t173; t123 * MDP(18) + t91 * MDP(19) - t92 * MDP(20) + (t144 * t175 + t82) * MDP(26) + (-t171 + (-t89 - t175) * t140) * MDP(27) + t151 + t159; t179; t149; 0.2e1 * t152 + t160; t173; t82 * MDP(26) - t83 * MDP(27) + t159; t155; t163; MDP(25) + t152; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
