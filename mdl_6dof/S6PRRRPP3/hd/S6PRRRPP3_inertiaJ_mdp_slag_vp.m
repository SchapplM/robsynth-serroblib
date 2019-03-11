% Calculate joint inertia matrix for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRRPP3_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:58:32
% EndTime: 2019-03-08 22:58:34
% DurationCPUTime: 0.64s
% Computational Cost: add. (575->183), mult. (1134->238), div. (0->0), fcn. (1081->8), ass. (0->64)
t158 = MDP(19) + MDP(23);
t172 = MDP(22) * pkin(9);
t179 = t158 + t172;
t129 = sin(qJ(3));
t178 = 0.2e1 * t129;
t157 = MDP(20) - MDP(25);
t176 = MDP(17) - t157;
t155 = MDP(21) + MDP(24);
t149 = MDP(18) - t155;
t154 = MDP(22) + MDP(26);
t126 = pkin(4) + qJ(6);
t131 = cos(qJ(4));
t128 = sin(qJ(4));
t169 = t128 * qJ(5);
t175 = -t126 * t131 - t169;
t174 = pkin(5) + pkin(9);
t132 = cos(qJ(3));
t173 = pkin(8) * t132;
t170 = sin(pkin(6));
t119 = qJ(5) * t131;
t168 = t128 * t129;
t167 = t129 * t131;
t166 = t132 * qJ(5);
t165 = pkin(4) * t168 + t129 * pkin(8);
t111 = -pkin(3) * t132 - pkin(9) * t129 - pkin(2);
t164 = -t131 * t111 + t128 * t173;
t102 = t128 * t111 + t131 * t173;
t122 = t128 ^ 2;
t124 = t131 ^ 2;
t163 = t122 + t124;
t162 = MDP(11) * t129;
t161 = MDP(13) * t131;
t107 = -pkin(3) + t175;
t160 = t107 * MDP(26);
t145 = -pkin(4) * t131 - t169;
t110 = -pkin(3) + t145;
t159 = t110 * MDP(22);
t156 = MDP(21) - MDP(18);
t121 = t132 * pkin(4);
t100 = t121 + t164;
t150 = sin(qJ(2)) * t170;
t148 = -pkin(4) * MDP(22) + MDP(20);
t146 = t170 * cos(qJ(2));
t99 = -t102 + t166;
t143 = t131 * MDP(14) - t128 * MDP(15);
t142 = -t164 * MDP(17) - t102 * MDP(18);
t141 = -t128 * MDP(24) - t131 * MDP(25);
t138 = -MDP(26) * t126 - MDP(25) + t148;
t113 = t174 * t131;
t137 = MDP(17) * pkin(3) + MDP(20) * t110 + MDP(23) * t113 - MDP(25) * t107;
t112 = t174 * t128;
t136 = -MDP(18) * pkin(3) - MDP(21) * t110 + MDP(23) * t112 - MDP(24) * t107;
t135 = -t128 * MDP(14) - t131 * MDP(15) - t113 * MDP(24) + t112 * MDP(25);
t133 = qJ(5) ^ 2;
t125 = cos(pkin(6));
t106 = t125 * t129 + t132 * t150;
t105 = -t125 * t132 + t129 * t150;
t103 = -qJ(5) * t167 + t165;
t98 = (qJ(6) * t128 - t119) * t129 + t165;
t96 = t106 * t131 - t128 * t146;
t95 = t106 * t128 + t131 * t146;
t92 = -pkin(5) * t168 - t99;
t90 = pkin(5) * t167 + qJ(6) * t132 + t100;
t1 = [MDP(1) + t154 * (t105 ^ 2 + t95 ^ 2 + t96 ^ 2); -MDP(4) * t150 + (t100 * t95 + t103 * t105 - t96 * t99) * MDP(22) + (t105 * t98 + t90 * t95 + t92 * t96) * MDP(26) + t158 * (t95 * t167 - t96 * t168) + t176 * (t105 * t168 + t132 * t95) - t149 * (-t105 * t167 - t132 * t96) + (MDP(10) * t132 + MDP(3) - t162) * t146; MDP(2) - 0.2e1 * pkin(2) * t162 + (t100 ^ 2 + t103 ^ 2 + t99 ^ 2) * MDP(22) + (t90 ^ 2 + t92 ^ 2 + t98 ^ 2) * MDP(26) + (0.2e1 * pkin(2) * MDP(10) + t132 * MDP(16) + (MDP(6) - t143) * t178) * t132 + 0.2e1 * (-MDP(20) * t100 + MDP(21) * t99 - MDP(24) * t92 + MDP(25) * t90 - t142) * t132 + ((t100 * t131 + t128 * t99) * MDP(19) + (-t128 * t92 + t131 * t90) * MDP(23) + (-MDP(24) * t131 + MDP(25) * t128) * t98 + (-t128 * MDP(20) - t131 * MDP(21)) * t103) * t178 + (MDP(12) * t124 - 0.2e1 * t128 * t161 + MDP(5) + 0.2e1 * (t128 * MDP(17) + t131 * MDP(18)) * pkin(8)) * t129 ^ 2; -t106 * MDP(11) + (t112 * t95 + t113 * t96) * MDP(26) + (t149 * t128 - t176 * t131 - MDP(10) + t159 + t160) * t105 + t179 * (t95 * t128 + t131 * t96); (t90 * t128 + t92 * t131) * MDP(23) + (t112 * t90 + t113 * t92) * MDP(26) + (t141 + t160) * t98 + (t131 * MDP(20) - t128 * MDP(21) + t159) * t103 + (-pkin(8) * MDP(11) + MDP(8) + (-t156 * t131 + (MDP(17) - MDP(20)) * t128) * pkin(9) + t135) * t132 + (MDP(7) - pkin(8) * MDP(10) + (-t122 + t124) * MDP(13) + (-MDP(17) * pkin(8) + t136) * t131 + (MDP(12) * t131 + MDP(18) * pkin(8) - t137) * t128) * t129 + (MDP(19) + t172) * (t100 * t128 - t131 * t99); MDP(9) + t122 * MDP(12) + (t163 * pkin(9) ^ 2 + t110 ^ 2) * MDP(22) + (t107 ^ 2 + t112 ^ 2 + t113 ^ 2) * MDP(26) + 0.2e1 * t163 * MDP(19) * pkin(9) + 0.2e1 * t137 * t131 + 0.2e1 * (t136 + t161) * t128; (t154 * qJ(5) - t149) * t96 + (-MDP(17) + t138) * t95; (0.2e1 * t121 + t164) * MDP(20) + (-pkin(4) * t100 - qJ(5) * t99) * MDP(22) - t100 * MDP(25) + (qJ(5) * t92 - t126 * t90) * MDP(26) + (-MDP(16) + (-qJ(6) - t126) * MDP(25)) * t132 + (t145 * MDP(19) + t175 * MDP(23) + t141 * pkin(5) + t143) * t129 + t142 + t155 * (-0.2e1 * t166 + t102); (-pkin(4) * t128 + t119) * MDP(19) + (-t126 * t128 + t119) * MDP(23) + (qJ(5) * t113 - t112 * t126) * MDP(26) + ((MDP(22) * qJ(5) + t156) * t131 + (-MDP(17) + t148) * t128) * pkin(9) - t135; MDP(16) - 0.2e1 * pkin(4) * MDP(20) + (pkin(4) ^ 2 + t133) * MDP(22) + 0.2e1 * t126 * MDP(25) + (t126 ^ 2 + t133) * MDP(26) + 0.2e1 * t155 * qJ(5); t154 * t95; MDP(22) * t100 + MDP(26) * t90 - t157 * t132 + t158 * t167; t112 * MDP(26) + t179 * t128; t138; t154; t96 * MDP(26); -MDP(23) * t168 - t132 * MDP(24) + MDP(26) * t92; MDP(23) * t131 + MDP(26) * t113; MDP(26) * qJ(5) + MDP(24); 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
