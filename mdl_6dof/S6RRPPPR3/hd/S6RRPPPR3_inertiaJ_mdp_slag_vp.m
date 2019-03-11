% Calculate joint inertia matrix for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPPR3_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:15:40
% EndTime: 2019-03-09 08:15:41
% DurationCPUTime: 0.43s
% Computational Cost: add. (460->135), mult. (712->168), div. (0->0), fcn. (618->6), ass. (0->60)
t118 = sin(qJ(6));
t120 = cos(qJ(6));
t116 = cos(pkin(9));
t119 = sin(qJ(2));
t121 = cos(qJ(2));
t145 = pkin(8) * t121;
t115 = sin(pkin(9));
t97 = -t121 * pkin(2) - t119 * qJ(3) - pkin(1);
t90 = t121 * pkin(3) - t97;
t83 = t119 * pkin(4) + t121 * qJ(5) + t90;
t99 = (pkin(7) - qJ(4)) * t119;
t77 = -t115 * t99 + t116 * t83;
t75 = t119 * pkin(5) + t116 * t145 + t77;
t78 = t115 * t83 + t116 * t99;
t76 = t115 * t145 + t78;
t92 = t120 * t115 + t118 * t116;
t84 = t92 * t121;
t91 = t118 * t115 - t120 * t116;
t85 = t91 * t121;
t152 = (-t118 * t76 + t120 * t75) * MDP(28) - (t118 * t75 + t120 * t76) * MDP(29) + t85 * MDP(25) + t84 * MDP(26);
t102 = t115 ^ 2 + t116 ^ 2;
t149 = t102 * MDP(21);
t148 = t102 * MDP(22);
t106 = t121 * qJ(4);
t100 = t121 * pkin(7) - t106;
t147 = t100 ^ 2;
t117 = qJ(3) + pkin(4);
t103 = t116 * pkin(5) + t117;
t146 = -0.2e1 * t103;
t122 = -pkin(2) - pkin(3);
t110 = -qJ(5) + t122;
t144 = pkin(8) - t110;
t88 = t91 * MDP(28);
t143 = -t92 * MDP(29) - t88;
t74 = -t77 * t115 + t78 * t116;
t142 = t74 * MDP(22);
t141 = t85 * MDP(23);
t140 = MDP(18) + t148;
t139 = t115 * MDP(20);
t138 = t116 * MDP(19);
t137 = t117 * MDP(22);
t136 = -pkin(2) * MDP(14) - MDP(11);
t135 = t122 * MDP(18) + MDP(16);
t134 = t78 * t115 + t77 * t116;
t131 = -t84 * MDP(28) + t85 * MDP(29);
t130 = -t92 * MDP(28) + t91 * MDP(29);
t129 = t138 - t139;
t128 = -MDP(19) * t115 - MDP(20) * t116;
t127 = -MDP(17) + t128;
t126 = t129 + t137;
t95 = t144 * t115;
t96 = t144 * t116;
t125 = -t92 * MDP(25) + t91 * MDP(26) + (t118 * t96 + t120 * t95) * MDP(28) - (t118 * t95 - t120 * t96) * MDP(29);
t124 = pkin(7) ^ 2;
t123 = qJ(3) ^ 2;
t114 = t121 ^ 2;
t113 = t119 ^ 2;
t87 = t102 * t110;
t86 = -t106 + (-pkin(5) * t115 + pkin(7)) * t121;
t1 = [(t114 * t124 + t97 ^ 2) * MDP(14) + (t90 ^ 2 + t99 ^ 2 + t147) * MDP(18) + (t77 ^ 2 + t78 ^ 2 + t147) * MDP(22) + MDP(1) + (0.2e1 * t84 * MDP(24) + t141) * t85 + (t124 * MDP(14) + MDP(27) + MDP(4)) * t113 + 0.2e1 * t131 * t86 + 0.2e1 * (t113 + t114) * MDP(12) * pkin(7) + 0.2e1 * (-pkin(1) * MDP(10) - t97 * MDP(13) + t90 * MDP(15) - t99 * MDP(17) + t77 * MDP(19) - t78 * MDP(20) + t121 * MDP(5) + t152) * t119 + 0.2e1 * (-t97 * MDP(11) - t90 * MDP(16) + t134 * MDP(21) + pkin(1) * MDP(9) + t127 * t100) * t121; -t74 * MDP(21) - t92 * t141 + (-t92 * t84 + t85 * t91) * MDP(24) + (-t103 * t84 - t86 * t91) * MDP(28) + (t103 * t85 - t86 * t92) * MDP(29) + t135 * t99 + t110 * t142 + (qJ(3) * MDP(18) + MDP(15) + t126) * t100 + (MDP(7) + t128 * t117 + (MDP(12) - MDP(17)) * qJ(3)) * t121 + (-pkin(2) * MDP(12) - t122 * MDP(17) + t128 * t110 + MDP(6) + t125) * t119 + ((qJ(3) * MDP(14) - MDP(10) + MDP(13)) * t121 + (-MDP(9) + t136) * t119) * pkin(7); MDP(8) + 0.2e1 * pkin(2) * MDP(11) + (pkin(2) ^ 2 + t123) * MDP(14) + (t122 ^ 2 + t123) * MDP(18) + t88 * t146 + t110 ^ 2 * t148 + (t137 + 0.2e1 * t138 - 0.2e1 * t139) * t117 + 0.2e1 * (MDP(13) + MDP(15)) * qJ(3) + 0.2e1 * t122 * MDP(16) - 0.2e1 * t87 * MDP(21) + (MDP(23) * t92 - 0.2e1 * t91 * MDP(24) + MDP(29) * t146) * t92; t99 * MDP(18) + t142 + (MDP(14) * pkin(7) + MDP(12) + t127 + t130) * t119; t87 * MDP(22) + t135 + t136 - t149; MDP(14) + t140; t90 * MDP(18) + t134 * MDP(22) + (-MDP(16) + t149) * t121 + (MDP(15) + t129 + t143) * t119; 0; 0; t140; t100 * MDP(22) + t128 * t121 + t131; t126 + t143; 0; 0; MDP(22); t119 * MDP(27) + t152; t125; t130; t143; 0; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
