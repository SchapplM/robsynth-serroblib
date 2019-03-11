% Calculate joint inertia matrix for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRPRP9_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:29:29
% EndTime: 2019-03-09 03:29:31
% DurationCPUTime: 0.66s
% Computational Cost: add. (764->176), mult. (1380->250), div. (0->0), fcn. (1359->6), ass. (0->70)
t162 = MDP(17) * qJ(4);
t112 = sin(pkin(9));
t113 = cos(pkin(9));
t144 = t112 ^ 2 + t113 ^ 2;
t161 = t144 * t162;
t136 = MDP(23) + MDP(25);
t135 = MDP(24) - MDP(27);
t114 = sin(qJ(5));
t116 = cos(qJ(5));
t100 = t112 * t116 + t113 * t114;
t139 = t113 * MDP(14);
t140 = t112 * MDP(15);
t150 = pkin(3) * MDP(17);
t121 = -t139 + t140 - t150;
t107 = -pkin(4) * t113 - pkin(3);
t145 = t116 * t113;
t99 = t112 * t114 - t145;
t83 = pkin(5) * t99 - qJ(6) * t100 + t107;
t160 = t83 * MDP(28) + t100 * t135 + t136 * t99 + t121;
t117 = cos(qJ(3));
t159 = 0.2e1 * t117;
t158 = -2 * MDP(19);
t157 = 2 * MDP(24);
t156 = 2 * MDP(25);
t155 = 2 * MDP(26);
t154 = 2 * MDP(27);
t153 = (pkin(1) * MDP(6));
t115 = sin(qJ(3));
t152 = pkin(5) * t115;
t151 = pkin(8) + qJ(4);
t118 = -pkin(1) - pkin(7);
t147 = t112 * t118;
t101 = pkin(3) * t115 - qJ(4) * t117 + qJ(2);
t96 = t113 * t101;
t84 = -pkin(8) * t113 * t117 + t96 + (pkin(4) - t147) * t115;
t148 = t112 * t117;
t146 = t115 * t118;
t90 = t112 * t101 + t113 * t146;
t86 = -pkin(8) * t148 + t90;
t79 = t114 * t84 + t116 * t86;
t149 = qJ(6) * t115;
t110 = t115 ^ 2;
t111 = t117 ^ 2;
t143 = -t110 - t111;
t142 = MDP(17) * t115;
t141 = t100 * MDP(18);
t138 = t118 * MDP(17);
t137 = MDP(17) + MDP(28);
t134 = t151 * t112;
t133 = t114 * t86 - t116 * t84;
t132 = -MDP(28) * pkin(5) - MDP(25);
t97 = pkin(4) * t148 - t117 * t118;
t130 = t144 * MDP(16);
t129 = -MDP(23) + t132;
t125 = MDP(28) * qJ(6) - t135;
t92 = t100 * t117;
t94 = -t114 * t148 + t117 * t145;
t124 = t94 * MDP(20) - t92 * MDP(21);
t123 = t100 * MDP(20) - t99 * MDP(21);
t122 = -t112 * MDP(14) - t113 * MDP(15);
t102 = t151 * t113;
t93 = t99 * t115;
t91 = t100 * t115;
t89 = -t112 * t146 + t96;
t88 = t116 * t102 - t114 * t134;
t87 = t102 * t114 + t116 * t134;
t80 = pkin(5) * t92 - qJ(6) * t94 + t97;
t77 = t133 - t152;
t76 = t149 + t79;
t1 = [MDP(1) + t111 * MDP(7) + (t111 * t118 ^ 2 + t89 ^ 2 + t90 ^ 2) * MDP(17) + t110 * MDP(22) + (t76 ^ 2 + t77 ^ 2 + t80 ^ 2) * MDP(28) + (t94 * MDP(18) + t158 * t92) * t94 + ((-2 * MDP(4) + t153) * pkin(1)) + (MDP(13) * t159 + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + 0.2e1 * (qJ(2) * MDP(12) - t117 * MDP(8) + t124) * t115 + 0.2e1 * (-t111 * t147 + t115 * t89) * MDP(14) + 0.2e1 * (-t111 * t113 * t118 - t115 * t90) * MDP(15) + 0.2e1 * (-t115 * t133 + t92 * t97) * MDP(23) + (-t115 * t79 + t94 * t97) * t157 + (-t115 * t77 + t80 * t92) * t156 + (-t76 * t92 + t77 * t94) * t155 + (t115 * t76 - t80 * t94) * t154 + (-t112 * t90 - t113 * t89) * MDP(16) * t159; MDP(4) - t153 + t111 * t138 + (t91 * t94 + t92 * t93) * MDP(26) + (-t117 * t80 - t76 * t93 + t77 * t91) * MDP(28) + (MDP(15) * t143 + t142 * t90) * t113 + (MDP(14) * t143 - t142 * t89) * t112 + t136 * (-t115 * t91 - t117 * t92) - t135 * (-t115 * t93 + t117 * t94); MDP(6) + (t110 * t144 + t111) * MDP(17) + (t91 ^ 2 + t93 ^ 2 + t111) * MDP(28); t94 * t141 + (-t100 * t92 - t94 * t99) * MDP(19) + (t107 * t92 + t97 * t99) * MDP(23) + (t100 * t97 + t107 * t94) * MDP(24) + (t80 * t99 + t83 * t92) * MDP(25) + (t100 * t77 - t76 * t99 + t87 * t94 - t88 * t92) * MDP(26) + (-t100 * t80 - t83 * t94) * MDP(27) + (t76 * t88 + t77 * t87 + t80 * t83) * MDP(28) + (MDP(9) + t122 * pkin(3) + (MDP(12) - t121) * t118) * t117 + (-t118 * MDP(13) + qJ(4) * t122 - t135 * t88 - t136 * t87 - MDP(10) + t123) * t115 + (MDP(16) + t162) * (-t112 * t89 + t90 * t113); (t100 * t91 + t93 * t99) * MDP(26) + (t87 * t91 - t88 * t93) * MDP(28) + (-MDP(13) + t130 + t161) * t115 + (MDP(12) - t160) * t117; MDP(11) + (t83 ^ 2 + t87 ^ 2 + t88 ^ 2) * MDP(28) + (0.2e1 * t139 - 0.2e1 * t140 + t150) * pkin(3) + 0.2e1 * (t107 * MDP(23) + t83 * MDP(25) - t88 * MDP(26)) * t99 + (-0.2e1 * t83 * MDP(27) + t107 * t157 + t155 * t87 + t158 * t99 + t141) * t100 + (0.2e1 * t130 + t161) * qJ(4); t80 * MDP(28) + t135 * t94 + t136 * t92 + (-t122 - t138) * t117; -t137 * t117; t160; t137; t115 * MDP(22) - t133 * MDP(23) - t79 * MDP(24) + (-t133 + 0.2e1 * t152) * MDP(25) + (-pkin(5) * t94 - qJ(6) * t92) * MDP(26) + (0.2e1 * t149 + t79) * MDP(27) + (-pkin(5) * t77 + qJ(6) * t76) * MDP(28) + t124; -t125 * t93 + t129 * t91; (-pkin(5) * t100 - qJ(6) * t99) * MDP(26) + t125 * t88 + t129 * t87 + t123; 0; MDP(22) + pkin(5) * t156 + qJ(6) * t154 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(28); -t115 * MDP(25) + t94 * MDP(26) + t77 * MDP(28); t91 * MDP(28); t100 * MDP(26) + t87 * MDP(28); 0; t132; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
