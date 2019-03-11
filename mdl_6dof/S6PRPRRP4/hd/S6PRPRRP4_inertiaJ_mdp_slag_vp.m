% Calculate joint inertia matrix for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRRP4_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:12:20
% EndTime: 2019-03-08 20:12:22
% DurationCPUTime: 0.45s
% Computational Cost: add. (689->144), mult. (1401->209), div. (0->0), fcn. (1571->10), ass. (0->69)
t106 = sin(pkin(11));
t108 = cos(pkin(11));
t111 = sin(qJ(4));
t152 = cos(qJ(4));
t96 = t152 * t106 + t111 * t108;
t159 = 0.2e1 * t96;
t158 = MDP(8) * qJ(3);
t110 = sin(qJ(5));
t104 = t110 ^ 2;
t113 = cos(qJ(5));
t105 = t113 ^ 2;
t137 = t104 + t105;
t157 = t137 * MDP(24);
t130 = MDP(26) * pkin(9) + MDP(24);
t135 = MDP(21) + MDP(23);
t134 = MDP(22) - MDP(25);
t156 = -t134 * t110 + t135 * t113 + MDP(14);
t101 = -pkin(3) * t108 - pkin(2);
t155 = 0.2e1 * t101;
t95 = t106 * t111 - t152 * t108;
t154 = pkin(5) * t95;
t153 = pkin(9) * t95;
t151 = pkin(2) * MDP(8);
t150 = pkin(8) + qJ(3);
t88 = pkin(4) * t95 - pkin(9) * t96 + t101;
t98 = t150 * t106;
t99 = t150 * t108;
t91 = -t111 * t98 + t152 * t99;
t79 = t110 * t88 + t113 * t91;
t148 = qJ(6) * t95;
t107 = sin(pkin(6));
t114 = cos(qJ(2));
t143 = t107 * t114;
t109 = cos(pkin(6));
t112 = sin(qJ(2));
t144 = t107 * t112;
t92 = -t106 * t144 + t108 * t109;
t93 = t106 * t109 + t108 * t144;
t85 = t111 * t92 + t152 * t93;
t83 = -t110 * t143 + t113 * t85;
t147 = t110 * t83;
t146 = t113 * t96;
t145 = t106 * MDP(6);
t142 = t108 * MDP(5);
t141 = t95 * MDP(20);
t140 = t96 * MDP(15);
t126 = -pkin(5) * t113 - qJ(6) * t110;
t97 = -pkin(4) + t126;
t139 = t97 * MDP(26);
t136 = MDP(17) * t113;
t132 = t110 * t91 - t113 * t88;
t131 = -MDP(26) * pkin(5) - MDP(23);
t129 = -pkin(4) * t96 - t153;
t128 = -t96 * t97 + t153;
t127 = MDP(21) - t131;
t125 = pkin(5) * t110 - qJ(6) * t113;
t122 = MDP(26) * qJ(6) - t134;
t121 = -t132 * MDP(21) - t79 * MDP(22);
t90 = t111 * t99 + t152 * t98;
t80 = t125 * t96 + t90;
t120 = -t90 * MDP(21) - t80 * MDP(23);
t119 = t90 * MDP(22) - t80 * MDP(25);
t118 = t113 * MDP(18) - t110 * MDP(19);
t117 = -t142 + t145 - t151;
t84 = t111 * t93 - t152 * t92;
t82 = t110 * t85 + t113 * t143;
t77 = t132 - t154;
t76 = t148 + t79;
t1 = [MDP(1) + (t107 ^ 2 * t114 ^ 2 + t92 ^ 2 + t93 ^ 2) * MDP(8) + (t82 ^ 2 + t83 ^ 2 + t84 ^ 2) * MDP(26); (t77 * t82 + t80 * t84) * MDP(26) + (t76 * MDP(26) - t134 * t95) * t83 + (-MDP(24) * t147 + (t82 * MDP(24) + t134 * t84) * t113) * t96 + (-t112 * MDP(4) + (-t95 * MDP(14) + MDP(3) - t117 - t140) * t114) * t107 + t135 * (t84 * t110 * t96 - t82 * t95) + (MDP(7) + t158) * (-t106 * t92 + t108 * t93); MDP(2) + t140 * t155 + (t76 ^ 2 + t77 ^ 2 + t80 ^ 2) * MDP(26) + (0.2e1 * t142 - 0.2e1 * t145 + t151) * pkin(2) + (t105 * MDP(16) - 0.2e1 * t110 * t136 + MDP(9)) * t96 ^ 2 + (MDP(14) * t155 + t141 + (-MDP(10) + t118) * t159) * t95 + 0.2e1 * (-t77 * MDP(23) + t76 * MDP(25) + t121) * t95 + ((t77 * MDP(24) + t119) * t113 + (-t76 * MDP(24) - t120) * t110) * t159 + (0.2e1 * MDP(7) + t158) * (t106 ^ 2 + t108 ^ 2) * qJ(3); -MDP(8) * t143 + (-t113 * t82 + t147) * MDP(26); (t110 * t76 - t113 * t77) * MDP(26) + (MDP(15) - t157) * t96 + t156 * t95 + t117; t137 * MDP(26) + MDP(8); -t85 * MDP(15) + (t139 - t156) * t84 + t130 * (t110 * t82 + t113 * t83); t80 * t139 - t95 * MDP(12) - t90 * MDP(14) - t91 * MDP(15) + (MDP(11) + (-t104 + t105) * MDP(17)) * t96 + (t95 * MDP(19) + t129 * MDP(22) + t128 * MDP(25) + t130 * t76 + t120) * t113 + (MDP(16) * t146 + t95 * MDP(18) + t129 * MDP(21) - t128 * MDP(23) + t130 * t77 + t119) * t110; 0; MDP(13) + t104 * MDP(16) + (t137 * pkin(9) ^ 2 + t97 ^ 2) * MDP(26) + 0.2e1 * pkin(9) * t157 + 0.2e1 * (MDP(21) * pkin(4) - MDP(23) * t97) * t113 + 0.2e1 * (-MDP(22) * pkin(4) - MDP(25) * t97 + t136) * t110; t122 * t83 - t127 * t82; t141 + (-t132 + 0.2e1 * t154) * MDP(23) + (0.2e1 * t148 + t79) * MDP(25) + (-pkin(5) * t77 + qJ(6) * t76) * MDP(26) + (t126 * MDP(24) + t118) * t96 + t121; t122 * t110 + t127 * t113; t110 * MDP(18) + t113 * MDP(19) - t125 * MDP(24) + (-t127 * t110 + t122 * t113) * pkin(9); MDP(20) + 0.2e1 * pkin(5) * MDP(23) + 0.2e1 * qJ(6) * MDP(25) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(26); t82 * MDP(26); -t95 * MDP(23) + MDP(24) * t146 + MDP(26) * t77; -t113 * MDP(26); t130 * t110; t131; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
