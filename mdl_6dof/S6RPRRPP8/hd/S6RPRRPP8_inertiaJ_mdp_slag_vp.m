% Calculate joint inertia matrix for
% S6RPRRPP8
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
%   see S6RPRRPP8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRRPP8_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:55:54
% EndTime: 2019-03-09 04:55:55
% DurationCPUTime: 0.57s
% Computational Cost: add. (534->171), mult. (867->213), div. (0->0), fcn. (684->4), ass. (0->62)
t137 = MDP(22) - MDP(27);
t136 = MDP(23) + MDP(26);
t107 = pkin(4) + qJ(6);
t111 = cos(qJ(4));
t109 = sin(qJ(4));
t148 = qJ(5) * t109;
t153 = -t107 * t111 - t148;
t152 = -(2 * pkin(4) * MDP(22)) + 0.2e1 * t107 * MDP(27) + MDP(18);
t112 = cos(qJ(3));
t151 = 0.2e1 * t112;
t150 = pkin(5) + pkin(8);
t149 = (pkin(1) * MDP(6));
t110 = sin(qJ(3));
t113 = -pkin(1) - pkin(7);
t145 = t110 * t113;
t93 = pkin(3) * t110 - pkin(8) * t112 + qJ(2);
t83 = -t109 * t145 + t111 * t93;
t84 = t109 * t93 + t111 * t145;
t102 = qJ(5) * t111;
t147 = t109 * t112;
t146 = t110 * qJ(5);
t144 = t111 * t112;
t125 = -pkin(4) * t111 - t148;
t94 = -pkin(3) + t125;
t143 = t94 * MDP(24);
t103 = t109 ^ 2;
t105 = t111 ^ 2;
t142 = t103 + t105;
t140 = t111 * MDP(15);
t139 = -MDP(20) + MDP(23);
t138 = MDP(21) + MDP(25);
t135 = MDP(24) + MDP(28);
t133 = t142 * pkin(8);
t132 = MDP(19) - t137;
t131 = -MDP(20) + t136;
t129 = MDP(24) * pkin(8) + MDP(21);
t128 = -pkin(4) * MDP(24) + MDP(22);
t127 = -t113 - t102;
t126 = t135 * t110;
t81 = -t146 - t84;
t124 = t131 * t111;
t123 = t83 * MDP(19) - t84 * MDP(20);
t122 = t111 * MDP(16) - t109 * MDP(17);
t121 = -t107 * MDP(28) - MDP(27) + t128;
t79 = -pkin(5) * t147 - t81;
t100 = pkin(4) * t147;
t80 = t100 + (qJ(6) * t109 + t127) * t112;
t85 = t112 * t127 + t100;
t120 = t85 * MDP(22) + t79 * MDP(25) - t80 * MDP(27);
t78 = pkin(5) * t144 - t107 * t110 - t83;
t119 = -t85 * MDP(23) + t78 * MDP(25) - t80 * MDP(26);
t87 = -pkin(3) + t153;
t96 = t150 * t111;
t118 = pkin(3) * MDP(19) + t94 * MDP(22) + t96 * MDP(25) - t87 * MDP(27);
t95 = t150 * t109;
t117 = -pkin(3) * MDP(20) - t94 * MDP(23) + t95 * MDP(25) - t87 * MDP(26);
t116 = t109 * MDP(16) + t111 * MDP(17) + t96 * MDP(26) - t95 * MDP(27);
t114 = qJ(5) ^ 2;
t106 = t112 ^ 2;
t104 = t110 ^ 2;
t82 = -pkin(4) * t110 - t83;
t1 = [MDP(1) + t104 * MDP(18) + (t81 ^ 2 + t82 ^ 2 + t85 ^ 2) * MDP(24) + (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) * MDP(28) + ((-2 * MDP(4) + t149) * pkin(1)) + (MDP(13) * t151 + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + ((t82 * MDP(21) + t119) * t111 + (t81 * MDP(21) - t120) * t109) * t151 + (t105 * MDP(14) - 0.2e1 * t109 * t140 + MDP(7) + 0.2e1 * (-t109 * MDP(19) - t111 * MDP(20)) * t113) * t106 + 0.2e1 * (qJ(2) * MDP(12) + (-MDP(8) + t122) * t112 + t82 * MDP(22) - t81 * MDP(23) + t79 * MDP(26) - t78 * MDP(27) + t123) * t110; -t149 + MDP(4) + (-t85 * MDP(24) - t80 * MDP(28)) * t112 + ((t109 * t82 - t111 * t81) * MDP(24) + (t109 * t78 + t111 * t79) * MDP(28)) * t110 + (-t132 * t109 + t124) * (t104 + t106); MDP(6) + t135 * (t104 * t142 + t106); t85 * t143 + (t78 * t95 + t79 * t96 + t80 * t87) * MDP(28) + (-t129 * t81 + t120) * t111 + (t129 * t82 + t119) * t109 + (-t113 * MDP(13) - MDP(10) + (t139 * t111 + (-MDP(19) + MDP(22)) * t109) * pkin(8) + t116) * t110 + (MDP(9) + t113 * MDP(12) + (-t103 + t105) * MDP(15) + (t113 * MDP(19) + t117) * t111 + (t111 * MDP(14) - t113 * MDP(20) - t118) * t109) * t112; (-t87 * MDP(28) + t109 * t131 + t111 * t132 + MDP(12) - t143) * t112 + (t138 * t142 - MDP(13) + (t109 * t95 + t111 * t96) * MDP(28) + MDP(24) * t133) * t110; MDP(11) + t103 * MDP(14) + (pkin(8) ^ 2 * t142 + t94 ^ 2) * MDP(24) + (t87 ^ 2 + t95 ^ 2 + t96 ^ 2) * MDP(28) + 0.2e1 * MDP(21) * t133 + 0.2e1 * t118 * t111 + 0.2e1 * (t117 + t140) * t109; (-pkin(4) * t82 - qJ(5) * t81) * MDP(24) + (qJ(5) * t79 - t107 * t78) * MDP(28) + t152 * t110 + (t125 * MDP(21) + t153 * MDP(25) + (-t109 * MDP(26) - t111 * MDP(27)) * pkin(5) + t122) * t112 + t123 - t137 * t83 + t136 * (0.2e1 * t146 + t84); t126 * t102 + (t124 + (-MDP(19) + t121) * t109) * t110; (-pkin(4) * t109 + t102) * MDP(21) + (-t107 * t109 + t102) * MDP(25) + (qJ(5) * t96 - t107 * t95) * MDP(28) + ((qJ(5) * MDP(24) + t139) * t111 + (-MDP(19) + t128) * t109) * pkin(8) + t116; ((pkin(4) ^ 2) + t114) * MDP(24) + (t107 ^ 2 + t114) * MDP(28) + 0.2e1 * t136 * qJ(5) + t152; t82 * MDP(24) + t78 * MDP(28) + t110 * t137 + t138 * t144; t109 * t126; MDP(28) * t95 + (MDP(25) + t129) * t109; t121; t135; -MDP(25) * t147 + t110 * MDP(26) + t79 * MDP(28); t111 * t110 * MDP(28); t111 * MDP(25) + MDP(28) * t96; MDP(28) * qJ(5) + MDP(26); 0; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
