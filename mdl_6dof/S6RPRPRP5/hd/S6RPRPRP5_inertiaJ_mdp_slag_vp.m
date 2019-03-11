% Calculate joint inertia matrix for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRP5_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:16:22
% EndTime: 2019-03-09 03:16:23
% DurationCPUTime: 0.65s
% Computational Cost: add. (1266->170), mult. (2359->240), div. (0->0), fcn. (2660->8), ass. (0->72)
t127 = sin(pkin(10));
t129 = cos(pkin(10));
t131 = sin(qJ(5));
t170 = cos(qJ(5));
t176 = -t127 * t131 + t170 * t129;
t154 = t127 ^ 2 + t129 ^ 2;
t175 = t154 * MDP(17);
t114 = t170 * t127 + t131 * t129;
t148 = MDP(25) - MDP(28);
t149 = MDP(24) + MDP(26);
t150 = t129 * MDP(15);
t151 = t127 * MDP(16);
t174 = -t148 * t114 + t149 * t176 + t150 - t151;
t130 = cos(pkin(9));
t122 = -pkin(2) * t130 - pkin(1);
t173 = 0.2e1 * t122;
t172 = 2 * MDP(27);
t171 = cos(qJ(3));
t169 = pkin(1) * MDP(7);
t128 = sin(pkin(9));
t132 = sin(qJ(3));
t113 = t128 * t132 - t171 * t130;
t168 = pkin(5) * t113;
t167 = pkin(7) + qJ(2);
t166 = pkin(8) + qJ(4);
t115 = t171 * t128 + t132 * t130;
t160 = t115 * t129;
t101 = pkin(3) * t113 - qJ(4) * t115 + t122;
t116 = t167 * t128;
t118 = t167 * t130;
t106 = -t132 * t116 + t171 * t118;
t92 = t129 * t101 - t106 * t127;
t89 = pkin(4) * t113 - pkin(8) * t160 + t92;
t161 = t115 * t127;
t93 = t127 * t101 + t129 * t106;
t91 = -pkin(8) * t161 + t93;
t85 = t131 * t89 + t170 * t91;
t165 = pkin(3) * MDP(18);
t97 = t176 * t115;
t164 = t176 * t97;
t163 = MDP(19) * t97;
t162 = qJ(6) * t113;
t158 = t128 * MDP(5);
t157 = t130 * MDP(4);
t96 = t114 * t115;
t156 = t96 * MDP(22);
t155 = t97 * MDP(21);
t152 = t113 * MDP(23);
t121 = -pkin(4) * t129 - pkin(3);
t146 = t166 * t127;
t145 = t131 * t91 - t170 * t89;
t144 = -MDP(29) * pkin(5) - MDP(26);
t143 = t154 * MDP(18);
t104 = t171 * t116 + t132 * t118;
t142 = -MDP(24) + t144;
t141 = t127 * t93 + t129 * t92;
t140 = -t127 * t92 + t129 * t93;
t139 = MDP(29) * qJ(6) - t148;
t138 = -t145 * MDP(24) - t85 * MDP(25);
t136 = t127 * MDP(15) + t129 * MDP(16);
t135 = t114 * MDP(21) + MDP(22) * t176;
t95 = pkin(4) * t161 + t104;
t117 = t166 * t129;
t109 = t114 ^ 2;
t105 = t170 * t117 - t131 * t146;
t103 = t117 * t131 + t170 * t146;
t100 = -pkin(5) * t176 - qJ(6) * t114 + t121;
t94 = t114 * t96;
t86 = t96 * pkin(5) - t97 * qJ(6) + t95;
t83 = t145 - t168;
t82 = t162 + t85;
t1 = [(t104 ^ 2 + t92 ^ 2 + t93 ^ 2) * MDP(18) + (t82 ^ 2 + t83 ^ 2 + t86 ^ 2) * MDP(29) + MDP(1) + (-0.2e1 * t96 * MDP(20) + t163) * t97 + (MDP(14) * t173 + MDP(8) * t115) * t115 + (0.2e1 * t157 - 0.2e1 * t158 + t169) * pkin(1) + (MDP(13) * t173 - 0.2e1 * t115 * MDP(9) + t152 + 0.2e1 * t155 - 0.2e1 * t156) * t113 + (-t82 * t96 + t83 * t97) * t172 + 0.2e1 * (t96 * MDP(24) + t97 * MDP(25)) * t95 + 0.2e1 * (t96 * MDP(26) - t97 * MDP(28)) * t86 + 0.2e1 * (t92 * MDP(15) - t93 * MDP(16) - t83 * MDP(26) + t82 * MDP(28) + t138) * t113 + 0.2e1 * (-t141 * MDP(17) + t136 * t104) * t115 + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t128 ^ 2 + t130 ^ 2) * qJ(2); -t157 + t158 - t169 + t141 * MDP(18) + (-t94 - t164) * MDP(27) + (t114 * t82 - t176 * t83) * MDP(29) + (MDP(14) - t175) * t115 + (MDP(13) + t174) * t113; MDP(7) + t143 + (t176 ^ 2 + t109) * MDP(29); t115 * MDP(10) - t104 * MDP(13) - t106 * MDP(14) + (-pkin(3) * t161 - t104 * t129) * MDP(15) + (-pkin(3) * t160 + t104 * t127) * MDP(16) + t140 * MDP(17) + (-pkin(3) * t104 + t140 * qJ(4)) * MDP(18) + t114 * t163 + (-t94 + t164) * MDP(20) + (t121 * t96 - t176 * t95) * MDP(24) + (t114 * t95 + t121 * t97) * MDP(25) + (t100 * t96 - t176 * t86) * MDP(26) + (t103 * t97 - t105 * t96 + t114 * t83 + t176 * t82) * MDP(27) + (-t100 * t97 - t114 * t86) * MDP(28) + (t100 * t86 + t103 * t83 + t105 * t82) * MDP(29) + (-t136 * qJ(4) - t149 * t103 - t148 * t105 - MDP(11) + t135) * t113; (-t103 * t176 + t105 * t114) * MDP(29); MDP(12) + t109 * MDP(19) + (t100 ^ 2 + t103 ^ 2 + t105 ^ 2) * MDP(29) + (0.2e1 * t150 - 0.2e1 * t151 + t165) * pkin(3) + (t103 * t114 + t105 * t176) * t172 + (t143 * qJ(4) + 0.2e1 * t175) * qJ(4) + 0.2e1 * (MDP(25) * t121 - MDP(28) * t100) * t114 - 0.2e1 * (-MDP(20) * t114 + MDP(24) * t121 + MDP(26) * t100) * t176; MDP(18) * t104 + MDP(29) * t86 + t136 * t115 + t148 * t97 + t149 * t96; 0; MDP(29) * t100 - t165 - t174; MDP(18) + MDP(29); t155 - t156 + t152 + (-t145 + 0.2e1 * t168) * MDP(26) + (-pkin(5) * t97 - qJ(6) * t96) * MDP(27) + (0.2e1 * t162 + t85) * MDP(28) + (-pkin(5) * t83 + qJ(6) * t82) * MDP(29) + t138; t139 * t114 - t142 * t176; (-pkin(5) * t114 + qJ(6) * t176) * MDP(27) + t139 * t105 + t142 * t103 + t135; 0; MDP(23) + 0.2e1 * pkin(5) * MDP(26) + 0.2e1 * qJ(6) * MDP(28) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(29); -t113 * MDP(26) + t97 * MDP(27) + MDP(29) * t83; -t176 * MDP(29); MDP(27) * t114 + MDP(29) * t103; 0; t144; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
