% Calculate joint inertia matrix for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPRP8_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:19:42
% EndTime: 2019-03-09 17:19:44
% DurationCPUTime: 0.89s
% Computational Cost: add. (772->193), mult. (1387->272), div. (0->0), fcn. (1295->6), ass. (0->76)
t138 = sin(qJ(5));
t140 = sin(qJ(2));
t139 = sin(qJ(3));
t141 = cos(qJ(5));
t169 = t141 * t139;
t142 = cos(qJ(3));
t170 = t140 * t142;
t108 = t138 * t170 - t140 * t169;
t113 = t138 * t139 + t141 * t142;
t109 = t113 * t140;
t143 = cos(qJ(2));
t133 = t143 * pkin(3);
t120 = -t143 * pkin(2) - t140 * pkin(8) - pkin(1);
t171 = t139 * t143;
t167 = pkin(7) * t171 - t142 * t120;
t104 = t133 + t167;
t96 = t143 * pkin(4) - pkin(9) * t170 + t104;
t174 = pkin(7) * t142;
t106 = t139 * t120 + t143 * t174;
t103 = -t143 * qJ(4) + t106;
t98 = t139 * t140 * pkin(9) + t103;
t155 = t138 * t98 - t141 * t96;
t90 = t138 * t96 + t141 * t98;
t186 = t109 * MDP(24) - t108 * MDP(25) - t155 * MDP(27) - t90 * MDP(28);
t185 = 2 * MDP(19);
t184 = -t167 * MDP(16) - t106 * MDP(17) - t186;
t183 = 0.2e1 * t143;
t181 = -t142 * pkin(3) - t139 * qJ(4);
t119 = -pkin(2) + t181;
t110 = t142 * pkin(4) - t119;
t180 = 0.2e1 * t110;
t179 = -2 * MDP(23);
t178 = 2 * MDP(29);
t177 = -pkin(3) - pkin(4);
t176 = pkin(8) - pkin(9);
t175 = pkin(3) * t139;
t173 = MDP(30) * pkin(5);
t172 = pkin(3) * MDP(21);
t168 = t142 * qJ(4);
t134 = t139 ^ 2;
t136 = t142 ^ 2;
t166 = t134 + t136;
t165 = qJ(4) * MDP(20);
t163 = t109 * MDP(22);
t117 = t138 * qJ(4) - t141 * t177;
t161 = t117 * MDP(27);
t118 = t141 * qJ(4) + t138 * t177;
t160 = t118 * MDP(28);
t159 = t138 * MDP(28);
t158 = t142 * MDP(12);
t157 = MDP(17) - MDP(20);
t156 = MDP(26) + MDP(15);
t121 = t176 * t139;
t122 = t176 * t142;
t101 = -t141 * t121 + t138 * t122;
t153 = t103 * t142 + t104 * t139;
t102 = t138 * t121 + t141 * t122;
t152 = t142 * MDP(13) - t139 * MDP(14);
t150 = t139 * MDP(18) - t142 * MDP(20);
t148 = t141 * MDP(27) + MDP(18) - t159;
t147 = -MDP(26) - t160 - t161;
t124 = t140 * t168;
t100 = t124 + (t177 * t139 - pkin(7)) * t140;
t114 = -t138 * t142 + t169;
t146 = t114 * MDP(24) - t113 * MDP(25) - t101 * MDP(27) - t102 * MDP(28);
t145 = -t139 * MDP(13) + t146;
t126 = pkin(8) * t171;
t116 = -pkin(5) - t117;
t107 = -t124 + (pkin(7) + t175) * t140;
t99 = t113 * pkin(5) + t110;
t93 = -t113 * qJ(6) + t102;
t92 = -t114 * qJ(6) - t101;
t91 = t108 * pkin(5) + t100;
t88 = -t108 * qJ(6) + t90;
t87 = t143 * pkin(5) - t109 * qJ(6) - t155;
t1 = [(t103 ^ 2 + t104 ^ 2 + t107 ^ 2) * MDP(21) + (t87 ^ 2 + t88 ^ 2 + t91 ^ 2) * MDP(30) + MDP(1) + t156 * t143 ^ 2 + (t108 * t179 + t163) * t109 + (-t88 * t108 - t87 * t109) * t178 + 0.2e1 * (t108 * MDP(27) + t109 * MDP(28)) * t100 + (t104 * MDP(18) - t103 * MDP(20) + pkin(1) * MDP(9) - t184) * t183 + (-0.2e1 * pkin(1) * MDP(10) + (MDP(5) - t152) * t183 + (-t103 * t139 + t104 * t142) * t185 + 0.2e1 * t150 * t107 + (t136 * MDP(11) - 0.2e1 * t139 * t158 + MDP(4) + 0.2e1 * (t139 * MDP(16) + t142 * MDP(17)) * pkin(7)) * t140) * t140; t126 * MDP(16) + (-t107 * t142 + t126) * MDP(18) + t153 * MDP(19) - t107 * t139 * MDP(20) + (t153 * pkin(8) + t107 * t119) * MDP(21) + t114 * t163 + (-t114 * t108 - t109 * t113) * MDP(23) + (t100 * t113 + t110 * t108) * MDP(27) + (t100 * t114 + t110 * t109) * MDP(28) + (-t93 * t108 - t92 * t109 - t88 * t113 - t87 * t114) * MDP(29) + (t87 * t92 + t88 * t93 + t91 * t99) * MDP(30) + (-pkin(7) * MDP(10) + MDP(7) + (t157 * pkin(8) - MDP(14)) * t142 + t145) * t143 + (MDP(6) - pkin(7) * MDP(9) + t142 * t139 * MDP(11) + (-t134 + t136) * MDP(12) + (-pkin(2) * t139 - t174) * MDP(16) + (-pkin(2) * t142 + pkin(7) * t139) * MDP(17) + t150 * t119) * t140; MDP(8) + t134 * MDP(11) + (t166 * pkin(8) ^ 2 + t119 ^ 2) * MDP(21) + t113 * MDP(27) * t180 + (t92 ^ 2 + t93 ^ 2 + t99 ^ 2) * MDP(30) + (MDP(22) * t114 + MDP(28) * t180 + t113 * t179) * t114 + (-t93 * t113 - t92 * t114) * t178 + t166 * pkin(8) * t185 + 0.2e1 * (pkin(2) * MDP(16) - t119 * MDP(18)) * t142 + 0.2e1 * (-pkin(2) * MDP(17) - t119 * MDP(20) + t158) * t139; (-0.2e1 * t133 - t167) * MDP(18) + t106 * MDP(20) + (-t104 * pkin(3) + t103 * qJ(4)) * MDP(21) + (-t118 * t108 - t116 * t109) * MDP(29) + (t87 * t116 + t88 * t118) * MDP(30) + (-MDP(15) + t147 - 0.2e1 * t165) * t143 + (t181 * MDP(19) + t152) * t140 + t184; t142 * MDP(14) + (t168 - t175) * MDP(19) + (-t118 * t113 - t116 * t114) * MDP(29) + (t92 * t116 + t93 * t118) * MDP(30) + ((MDP(21) * qJ(4) - t157) * t142 + (-MDP(16) - MDP(18) - t172) * t139) * pkin(8) - t145; 0.2e1 * pkin(3) * MDP(18) + 0.2e1 * t165 + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(21) + (t116 ^ 2 + t118 ^ 2) * MDP(30) + 0.2e1 * t161 + 0.2e1 * t160 + t156; MDP(19) * t170 + t104 * MDP(21) + (-t138 * t108 - t141 * t109) * MDP(29) + (t88 * t138 + t87 * t141) * MDP(30) + t148 * t143; (-t138 * t113 - t141 * t114) * MDP(29) + (t93 * t138 + t92 * t141) * MDP(30) + (pkin(8) * MDP(21) + MDP(19)) * t139; -t172 + (t116 * t141 + t118 * t138) * MDP(30) - t148; MDP(21) + (t138 ^ 2 + t141 ^ 2) * MDP(30); t143 * MDP(26) + (-t109 * MDP(29) + t87 * MDP(30)) * pkin(5) + t186; (-t114 * MDP(29) + t92 * MDP(30)) * pkin(5) + t146; t116 * t173 + t147; -t159 + (MDP(27) + t173) * t141; MDP(30) * pkin(5) ^ 2 + MDP(26); t91 * MDP(30); t99 * MDP(30); 0; 0; 0; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
