% Calculate joint inertia matrix for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP4_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:12:51
% EndTime: 2021-01-16 03:12:53
% DurationCPUTime: 0.49s
% Computational Cost: add. (437->155), mult. (826->202), div. (0->0), fcn. (776->8), ass. (0->70)
t154 = pkin(4) + pkin(8);
t108 = sin(qJ(5));
t152 = 0.2e1 * t108;
t107 = cos(pkin(6));
t109 = sin(qJ(3));
t112 = cos(qJ(3));
t106 = sin(pkin(6));
t110 = sin(qJ(2));
t144 = t106 * t110;
t88 = t107 * t109 + t112 * t144;
t153 = t88 ^ 2;
t114 = -pkin(3) - pkin(9);
t151 = t109 * pkin(5);
t150 = MDP(25) * pkin(5);
t149 = MDP(26) * pkin(5);
t148 = pkin(3) * MDP(15);
t147 = pkin(8) * MDP(15);
t95 = t154 * t109;
t146 = t108 * t95;
t87 = -t107 * t112 + t109 * t144;
t145 = t87 * t109;
t113 = cos(qJ(2));
t143 = t106 * t113;
t142 = t108 * t112;
t111 = cos(qJ(5));
t141 = t111 * t112;
t99 = t108 * pkin(5) + qJ(4);
t140 = t99 * MDP(23);
t139 = t99 * MDP(26);
t138 = -qJ(6) + t114;
t96 = t154 * t112;
t137 = qJ(4) * MDP(15);
t136 = t111 * MDP(17);
t135 = t111 * MDP(19);
t134 = t111 * MDP(24);
t133 = pkin(8) ^ 2 * MDP(15);
t132 = -MDP(11) + MDP(14);
t131 = MDP(21) + MDP(23);
t130 = MDP(22) + MDP(24);
t128 = -qJ(4) * t109 - pkin(2);
t90 = t112 * t114 + t128;
t79 = -t108 * t90 + t111 * t95;
t127 = MDP(12) + t147;
t126 = MDP(13) - t148;
t125 = MDP(23) + t149;
t124 = t130 * t108;
t123 = MDP(21) * t114 + MDP(18);
t122 = qJ(6) * t142 + t79;
t121 = -MDP(10) + t126;
t120 = MDP(21) + t125;
t78 = t146 + (-qJ(6) * t112 + t90) * t111;
t119 = t79 * MDP(21) - (t111 * t90 + t146) * MDP(22) - t78 * MDP(24);
t86 = pkin(5) * t141 + t96;
t118 = t96 * MDP(21) + t86 * MDP(23) - t78 * MDP(25);
t77 = t122 + t151;
t117 = t96 * MDP(22) + t86 * MDP(24) - t77 * MDP(25);
t92 = t138 * t108;
t116 = (-MDP(22) * t114 - MDP(19)) * t108 - t92 * MDP(24);
t105 = t112 ^ 2;
t104 = t111 ^ 2;
t103 = t109 ^ 2;
t102 = t108 ^ 2;
t97 = t102 + t104;
t94 = -t112 * pkin(3) + t128;
t93 = t138 * t111;
t83 = -t87 * t108 + t111 * t143;
t82 = t108 * t143 + t87 * t111;
t81 = t92 * t108 + t93 * t111;
t74 = -t83 * t108 + t82 * t111;
t1 = [MDP(1) + (t106 ^ 2 * t113 ^ 2 + t87 ^ 2 + t153) * MDP(15) + (t82 ^ 2 + t83 ^ 2 + t153) * MDP(26); MDP(12) * t145 + (t82 * t77 - t83 * t78 + t88 * t86) * MDP(26) + t130 * (t83 * t109 - t142 * t88) + t131 * (t82 * t109 + t141 * t88) + (t88 * MDP(12) + (t108 * t82 + t111 * t83) * MDP(25)) * t112 + (t88 * t112 + t145) * t147 + (-t110 * MDP(4) + (-MDP(15) * t94 + MDP(3) + (MDP(10) - MDP(13)) * t112 + t132 * t109) * t113) * t106; MDP(2) + t94 ^ 2 * MDP(15) + (t77 ^ 2 + t78 ^ 2 + t86 ^ 2) * MDP(26) + (t102 * MDP(16) + t136 * t152 + t133) * t105 + (MDP(20) + MDP(5) + t133) * t103 + 0.2e1 * (-pkin(2) * MDP(11) - t94 * MDP(14)) * t109 + 0.2e1 * (t77 * MDP(23) + t119) * t109 + 0.2e1 * (t103 + t105) * MDP(12) * pkin(8) + 0.2e1 * (pkin(2) * MDP(10) + t94 * MDP(13) + (-t108 * MDP(18) + MDP(6) - t135) * t109 - t108 * t117 + t111 * t118) * t112; -t74 * MDP(25) + (t82 * t93 - t83 * t92) * MDP(26) + t121 * t87 + (t108 * t131 + t111 * t130 + t132 + t137 + t139) * t88; (t77 * t93 + t78 * t92 + t86 * t99) * MDP(26) + t117 * t111 + t118 * t108 + (-pkin(3) * MDP(12) + t93 * MDP(23) + pkin(8) * t121 + t111 * t123 + MDP(7) + t116) * t109 + (MDP(8) + (t102 - t104) * MDP(17) + (-t92 * MDP(25) + t140) * t111 + (-t111 * MDP(16) - t99 * MDP(24) + t93 * MDP(25)) * t108 + t132 * pkin(8) + (t111 * MDP(21) - t108 * MDP(22) + t127) * qJ(4)) * t112; MDP(9) + t104 * MDP(16) + 0.2e1 * t99 * t134 - 0.2e1 * t81 * MDP(25) + (t92 ^ 2 + t93 ^ 2 + t99 ^ 2) * MDP(26) + (-t136 + t140) * t152 + (-0.2e1 * MDP(13) + t148) * pkin(3) + (MDP(21) * t152 + 0.2e1 * t111 * MDP(22) + 0.2e1 * MDP(14) + t137) * qJ(4); t87 * MDP(15) + t74 * MDP(26); (t78 * t108 + t77 * t111) * MDP(26) + (t111 * t131 - t124 + t127) * t109; -t97 * MDP(25) + t81 * MDP(26) + t126; t97 * MDP(26) + MDP(15); t120 * t82 + t130 * t83; t109 * MDP(20) + (t122 + 0.2e1 * t151) * MDP(23) + t77 * t149 + (-t135 + (-MDP(18) + t150) * t108) * t112 + t119; t125 * t93 + (t123 - t150) * t111 + t116; t111 * t120 - t124; MDP(20) + (0.2e1 * MDP(23) + t149) * pkin(5); t88 * MDP(26); t86 * MDP(26) + (t111 * MDP(23) - t108 * MDP(24)) * t112; t108 * MDP(23) + t134 + t139; 0; 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
