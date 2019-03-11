% Calculate joint inertia matrix for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRRPP7_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:52:14
% EndTime: 2019-03-09 04:52:16
% DurationCPUTime: 0.58s
% Computational Cost: add. (534->171), mult. (867->210), div. (0->0), fcn. (685->4), ass. (0->62)
t137 = MDP(21) + MDP(25);
t107 = sin(qJ(4));
t111 = pkin(4) + pkin(5);
t109 = cos(qJ(4));
t146 = t109 * qJ(5);
t153 = t111 * t107 - t146;
t152 = (2 * pkin(4) * MDP(21)) + 0.2e1 * t111 * MDP(25) + MDP(18);
t110 = cos(qJ(3));
t151 = 0.2e1 * t110;
t150 = (pkin(1) * MDP(6));
t149 = pkin(8) - qJ(6);
t108 = sin(qJ(3));
t113 = -pkin(1) - pkin(7);
t147 = t108 * t113;
t93 = t108 * pkin(3) - t110 * pkin(8) + qJ(2);
t148 = t107 * t147 - t109 * t93;
t84 = t107 * t93 + t109 * t147;
t145 = t109 * t110;
t80 = (t113 - t153) * t110;
t144 = t80 * MDP(28);
t131 = t107 * qJ(5) + pkin(3);
t87 = t111 * t109 + t131;
t143 = t87 * MDP(28);
t94 = -t109 * pkin(4) - t131;
t142 = t94 * MDP(24);
t103 = t107 ^ 2;
t105 = t109 ^ 2;
t141 = t103 + t105;
t139 = t109 * MDP(15);
t138 = -MDP(20) + MDP(23);
t136 = MDP(22) - MDP(27);
t135 = MDP(23) + MDP(26);
t134 = MDP(24) + MDP(28);
t102 = t108 * qJ(5);
t133 = 0.2e1 * t102 + t84;
t81 = t102 + t84;
t132 = t141 * pkin(8);
t130 = MDP(19) + t137;
t129 = -MDP(20) + t135;
t127 = -pkin(4) * MDP(24) - MDP(21);
t126 = MDP(24) * pkin(8) + MDP(22);
t125 = t134 * t108;
t124 = -pkin(4) * t107 + t146;
t123 = t129 * t109;
t122 = -t148 * MDP(19) - t84 * MDP(20);
t121 = -t111 * MDP(28) - MDP(25) + t127;
t99 = t107 * t110 * qJ(6);
t79 = t99 + t81;
t85 = (-t113 - t124) * t110;
t120 = -t85 * MDP(21) + t80 * MDP(25) - t79 * MDP(27);
t78 = -qJ(6) * t145 - t111 * t108 + t148;
t119 = -t85 * MDP(23) + t80 * MDP(26) - t78 * MDP(27);
t96 = t149 * t109;
t118 = pkin(3) * MDP(19) - t94 * MDP(21) + t87 * MDP(25) - t96 * MDP(27);
t95 = t149 * t107;
t117 = -pkin(3) * MDP(20) - t94 * MDP(23) + t87 * MDP(26) - t95 * MDP(27);
t116 = t107 * MDP(16) + t109 * MDP(17) - t95 * MDP(25) + t96 * MDP(26);
t114 = qJ(5) ^ 2;
t106 = t110 ^ 2;
t104 = t108 ^ 2;
t82 = -t108 * pkin(4) + t148;
t1 = [MDP(1) + t104 * MDP(18) + (t81 ^ 2 + t82 ^ 2 + t85 ^ 2) * MDP(24) + (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) * MDP(28) + ((-2 * MDP(4) + t150) * pkin(1)) + (MDP(13) * t151 + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + ((t82 * MDP(22) + t119) * t109 + (-t81 * MDP(22) - t120) * t107) * t151 + (t105 * MDP(14) - 0.2e1 * t107 * t139 + MDP(7) + 0.2e1 * (-t107 * MDP(19) - t109 * MDP(20)) * t113) * t106 + 0.2e1 * (qJ(2) * MDP(12) + (t109 * MDP(16) - t107 * MDP(17) - MDP(8)) * t110 - t82 * MDP(21) + t81 * MDP(23) - t78 * MDP(25) + t79 * MDP(26) + t122) * t108; -t150 + MDP(4) + (-t85 * MDP(24) + t144) * t110 + ((t82 * t107 + t81 * t109) * MDP(24) + (t78 * t107 + t79 * t109) * MDP(28)) * t108 + (-t130 * t107 + t123) * (t104 + t106); MDP(6) + t134 * (t141 * t104 + t106); t85 * t142 + (t78 * t95 + t79 * t96 + t80 * t87) * MDP(28) + (t126 * t81 + t120) * t109 + (t126 * t82 + t119) * t107 + (-t113 * MDP(13) - MDP(10) + (t138 * t109 + (-MDP(19) - MDP(21)) * t107) * pkin(8) + t116) * t108 + (MDP(9) + t113 * MDP(12) + (-t103 + t105) * MDP(15) + (t113 * MDP(19) + t117) * t109 + (t109 * MDP(14) - t113 * MDP(20) - t118) * t107) * t110; (t107 * t129 + t109 * t130 + MDP(12) - t142 + t143) * t110 + (t136 * t141 - MDP(13) + (t95 * t107 + t96 * t109) * MDP(28) + MDP(24) * t132) * t108; MDP(11) + t103 * MDP(14) + (t141 * pkin(8) ^ 2 + t94 ^ 2) * MDP(24) + (t87 ^ 2 + t95 ^ 2 + t96 ^ 2) * MDP(28) + 0.2e1 * MDP(22) * t132 + 0.2e1 * t118 * t109 + 0.2e1 * (t117 + t139) * t107; t133 * MDP(23) + (-t82 * pkin(4) + t81 * qJ(5)) * MDP(24) + (t99 + t133) * MDP(26) + (t79 * qJ(5) - t78 * t111) * MDP(28) + t152 * t108 + ((-t136 * qJ(5) - MDP(17)) * t107 + (-pkin(4) * MDP(22) + qJ(6) * MDP(25) + t111 * MDP(27) + MDP(16)) * t109) * t110 + t122 - t137 * t148; t125 * t146 + (t123 + (-MDP(19) + t121) * t107) * t108; t124 * MDP(22) + t153 * MDP(27) + (t96 * qJ(5) - t95 * t111) * MDP(28) + ((qJ(5) * MDP(24) + t138) * t109 + (-MDP(19) + t127) * t107) * pkin(8) + t116; ((pkin(4) ^ 2) + t114) * MDP(24) + (t111 ^ 2 + t114) * MDP(28) + 0.2e1 * t135 * qJ(5) + t152; t82 * MDP(24) + t78 * MDP(28) - t137 * t108 + t136 * t145; t107 * t125; t95 * MDP(28) + (-MDP(27) + t126) * t107; t121; t134; t144 + (-t107 * MDP(25) + t109 * MDP(26)) * t110; t110 * MDP(28); t109 * MDP(25) + t107 * MDP(26) + t143; 0; 0; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
