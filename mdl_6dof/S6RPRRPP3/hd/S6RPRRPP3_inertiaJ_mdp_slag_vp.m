% Calculate joint inertia matrix for
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRRPP3_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:37:00
% EndTime: 2019-03-09 04:37:02
% DurationCPUTime: 0.57s
% Computational Cost: add. (559->169), mult. (936->219), div. (0->0), fcn. (762->6), ass. (0->62)
t134 = MDP(21) + MDP(24);
t111 = pkin(4) + qJ(6);
t115 = cos(qJ(4));
t113 = sin(qJ(4));
t145 = qJ(5) * t113;
t151 = -t111 * t115 - t145;
t150 = pkin(5) + pkin(8);
t116 = cos(qJ(3));
t109 = sin(pkin(9));
t99 = pkin(1) * t109 + pkin(7);
t146 = t116 * t99;
t110 = cos(pkin(9));
t100 = -pkin(1) * t110 - pkin(2);
t114 = sin(qJ(3));
t89 = -pkin(3) * t116 - pkin(8) * t114 + t100;
t149 = t113 * t146 - t115 * t89;
t83 = t113 * t89 + t115 * t146;
t144 = t113 * t114;
t143 = t114 * t115;
t97 = qJ(5) * t143;
t148 = -pkin(4) * t144 + t97;
t147 = pkin(4) * MDP(22);
t142 = t116 * qJ(5);
t128 = -pkin(4) * t115 - t145;
t93 = -pkin(3) + t128;
t141 = t93 * MDP(22);
t105 = t113 ^ 2;
t107 = t115 ^ 2;
t140 = t105 + t107;
t139 = MDP(13) * t115;
t138 = MDP(19) + MDP(23);
t137 = -MDP(20) + MDP(25);
t136 = MDP(20) - MDP(17);
t135 = MDP(21) - MDP(18);
t133 = MDP(22) + MDP(26);
t84 = t114 * t99 - t148;
t104 = t116 * pkin(4);
t80 = t104 + t149;
t131 = t140 * pkin(8);
t130 = -MDP(18) + t134;
t129 = MDP(22) * pkin(8) + MDP(19);
t79 = t142 - t83;
t127 = -MDP(26) * t111 - t137;
t126 = -t149 * MDP(17) - t83 * MDP(18);
t125 = t115 * MDP(14) - t113 * MDP(15);
t78 = -pkin(5) * t144 - t79;
t81 = qJ(6) * t144 + t84;
t124 = t84 * MDP(20) + t78 * MDP(23) - t81 * MDP(25);
t77 = pkin(5) * t143 + qJ(6) * t116 + t80;
t123 = -t84 * MDP(21) + t77 * MDP(23) - t81 * MDP(24);
t88 = -pkin(3) + t151;
t96 = t150 * t115;
t122 = MDP(17) * pkin(3) + MDP(20) * t93 + MDP(23) * t96 - MDP(25) * t88;
t95 = t150 * t113;
t121 = -MDP(18) * pkin(3) - MDP(21) * t93 + MDP(23) * t95 - MDP(24) * t88;
t120 = -t113 * MDP(14) - t115 * MDP(15) - t96 * MDP(24) + t95 * MDP(25);
t117 = qJ(5) ^ 2;
t108 = t116 ^ 2;
t106 = t114 ^ 2;
t103 = qJ(5) * t115;
t101 = t107 * t114;
t1 = [MDP(1) - 0.2e1 * t100 * t116 * MDP(10) + t108 * MDP(16) + (t79 ^ 2 + t80 ^ 2 + t84 ^ 2) * MDP(22) + (t77 ^ 2 + t78 ^ 2 + t81 ^ 2) * MDP(26) + (t109 ^ 2 + t110 ^ 2) * MDP(4) * pkin(1) ^ 2 + 0.2e1 * (-t80 * MDP(20) + t79 * MDP(21) - t78 * MDP(24) + t77 * MDP(25) - t126) * t116 + (t107 * MDP(12) - 0.2e1 * t113 * t139 + MDP(5) + 0.2e1 * (t113 * MDP(17) + t115 * MDP(18)) * t99) * t106 + 0.2e1 * (t100 * MDP(11) + (MDP(6) - t125) * t116 + (t80 * MDP(19) + t123) * t115 + (t79 * MDP(19) - t124) * t113) * t114; (-MDP(22) * t84 - MDP(26) * t81) * t116 + ((t113 * t80 - t115 * t79) * MDP(22) + (t113 * t77 + t115 * t78) * MDP(26)) * t114; MDP(4) + t133 * (t140 * t106 + t108); t101 * MDP(13) + t84 * t141 + (t77 * t95 + t78 * t96 + t81 * t88) * MDP(26) + (-t129 * t79 + t124) * t115 + (t129 * t80 + t123) * t113 + (-t99 * MDP(11) + MDP(8) + (-t136 * t113 - t135 * t115) * pkin(8) + t120) * t116 + (-t99 * MDP(10) - t105 * MDP(13) + MDP(7) + (-MDP(17) * t99 + t121) * t115 + (MDP(12) * t115 + MDP(18) * t99 - t122) * t113) * t114; t138 * (t105 * t114 + t101) + (-MDP(11) + (t113 * t95 + t115 * t96) * MDP(26) + MDP(22) * t131) * t114 + (-t141 - t88 * MDP(26) + MDP(10) + (MDP(25) - t136) * t115 + t130 * t113) * t116; MDP(9) + t105 * MDP(12) + (t140 * pkin(8) ^ 2 + t93 ^ 2) * MDP(22) + (t88 ^ 2 + t95 ^ 2 + t96 ^ 2) * MDP(26) + 0.2e1 * MDP(19) * t131 + 0.2e1 * t122 * t115 + 0.2e1 * (t121 + t139) * t113; (0.2e1 * t104 + t149) * MDP(20) + (-pkin(4) * t80 - qJ(5) * t79) * MDP(22) - t80 * MDP(25) + (qJ(5) * t78 - t111 * t77) * MDP(26) + (-MDP(16) + (-qJ(6) - t111) * MDP(25)) * t116 + (t128 * MDP(19) + t151 * MDP(23) + (-t113 * MDP(24) - t115 * MDP(25)) * pkin(5) + t125) * t114 + t126 + t134 * (-0.2e1 * t142 + t83); t148 * MDP(22) + t97 * MDP(26) + (t130 * t115 + (-MDP(17) + t127) * t113) * t114; (-pkin(4) * t113 + t103) * MDP(19) + (-t111 * t113 + t103) * MDP(23) + (qJ(5) * t96 - t111 * t95) * MDP(26) + ((MDP(22) * qJ(5) + t135) * t115 + (t136 - t147) * t113) * pkin(8) - t120; MDP(16) - 0.2e1 * pkin(4) * MDP(20) + (pkin(4) ^ 2 + t117) * MDP(22) + 0.2e1 * t111 * MDP(25) + (t111 ^ 2 + t117) * MDP(26) + 0.2e1 * t134 * qJ(5); MDP(22) * t80 + MDP(26) * t77 + t137 * t116 + t138 * t143; t133 * t144; MDP(26) * t95 + (MDP(23) + t129) * t113; t127 - t147; t133; -MDP(23) * t144 - t116 * MDP(24) + MDP(26) * t78; MDP(26) * t143; MDP(23) * t115 + MDP(26) * t96; MDP(26) * qJ(5) + MDP(24); 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
