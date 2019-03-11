% Calculate joint inertia matrix for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PPRRRP2_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:58:10
% EndTime: 2019-03-08 18:58:12
% DurationCPUTime: 0.44s
% Computational Cost: add. (564->148), mult. (1394->231), div. (0->0), fcn. (1587->12), ass. (0->63)
t114 = sin(qJ(4));
t153 = 0.2e1 * t114;
t152 = MDP(23) * pkin(10) + MDP(21);
t151 = MDP(18) + MDP(20);
t113 = sin(qJ(5));
t150 = pkin(9) * t113;
t117 = cos(qJ(4));
t100 = -t117 * pkin(4) - t114 * pkin(10) - pkin(3);
t116 = cos(qJ(5));
t92 = t116 * t117 * pkin(9) + t113 * t100;
t108 = sin(pkin(7));
t115 = sin(qJ(3));
t148 = t108 * t115;
t118 = cos(qJ(3));
t147 = t108 * t118;
t110 = cos(pkin(12));
t111 = cos(pkin(7));
t146 = t110 * t111;
t145 = t113 * t114;
t144 = t113 * t117;
t143 = t114 * t116;
t104 = t113 ^ 2;
t106 = t116 ^ 2;
t142 = t104 + t106;
t141 = MDP(21) * t114;
t140 = t113 * MDP(15);
t139 = t114 * MDP(12);
t138 = t116 * MDP(14);
t137 = (MDP(19) - MDP(22));
t134 = -MDP(23) * pkin(5) - MDP(20);
t133 = -t117 * MDP(11) - MDP(4);
t132 = -2 * qJ(6) * MDP(22) - MDP(17);
t131 = -MDP(18) + t134;
t130 = -t116 * pkin(5) - t113 * qJ(6);
t129 = -pkin(5) * t113 + t116 * qJ(6);
t88 = -t117 * qJ(6) + t92;
t98 = t116 * t100;
t89 = -t98 + (pkin(5) + t150) * t117;
t126 = t89 * t113 + t88 * t116;
t99 = -pkin(4) + t130;
t125 = pkin(4) * MDP(18) - t99 * MDP(20);
t124 = -pkin(4) * MDP(19) - t99 * MDP(22);
t123 = MDP(23) * qJ(6) - t137;
t122 = (-pkin(9) * t144 + t98) * MDP(18) - t92 * MDP(19);
t121 = t116 * MDP(15) - t113 * MDP(16);
t120 = t99 * MDP(23) + t137 * t113 - t151 * t116 - MDP(11);
t112 = cos(pkin(6));
t109 = sin(pkin(6));
t107 = sin(pkin(12));
t101 = pkin(10) * t144;
t96 = t114 * t111 + t117 * t148;
t95 = -t117 * t111 + t114 * t148;
t94 = -t109 * t110 * t108 + t112 * t111;
t93 = (pkin(9) - t129) * t114;
t87 = -t113 * t147 + t116 * t96;
t86 = t113 * t96 + t116 * t147;
t85 = t112 * t148 + (t107 * t118 + t115 * t146) * t109;
t84 = -t112 * t147 + (t107 * t115 - t118 * t146) * t109;
t81 = t94 * t114 + t85 * t117;
t80 = t85 * t114 - t94 * t117;
t78 = t84 * t113 + t81 * t116;
t77 = t81 * t113 - t84 * t116;
t1 = [MDP(1) + (t112 ^ 2 + (t107 ^ 2 + t110 ^ 2) * t109 ^ 2) * MDP(2) + (t77 ^ 2 + t78 ^ 2 + t80 ^ 2) * MDP(23); t112 * MDP(2) + (t77 * t86 + t78 * t87 + t80 * t95) * MDP(23); MDP(2) + (t86 ^ 2 + t87 ^ 2 + t95 ^ 2) * MDP(23); -t85 * MDP(5) + (t77 * t89 + t78 * t88 + t80 * t93) * MDP(23) + t133 * t84 + t137 * (t78 * t117 + t80 * t143) + (t84 * MDP(12) + (-t113 * t78 + t116 * t77) * MDP(21)) * t114 + t151 * (t77 * t117 + t80 * t145); (t86 * t89 + t87 * t88 + t95 * t93) * MDP(23) + t137 * (t87 * t117 + t95 * t143) + (-t113 * t87 + t116 * t86) * t141 + (-t115 * MDP(5) + (-t133 - t139) * t118) * t108 + t151 * (t86 * t117 + t95 * t145); MDP(3) - 0.2e1 * pkin(3) * t139 + (t88 ^ 2 + t89 ^ 2 + t93 ^ 2) * MDP(23) + (0.2e1 * pkin(3) * MDP(11) + t117 * MDP(17) + (MDP(7) - t121) * t153) * t117 + 0.2e1 * (t89 * MDP(20) - t88 * MDP(22) - t122) * t117 + ((-t113 * t88 + t116 * t89) * MDP(21) + (t113 * MDP(20) - t116 * MDP(22)) * t93) * t153 + (t106 * MDP(13) - 0.2e1 * t113 * t138 + MDP(6) + 0.2e1 * (t113 * MDP(18) + t116 * MDP(19)) * pkin(9)) * t114 ^ 2; -t81 * MDP(12) + t120 * t80 + t152 * (t77 * t113 + t78 * t116); -t96 * MDP(12) + t120 * t95 + t152 * (t86 * t113 + t87 * t116); t101 * MDP(18) + (-t93 * t116 + t101) * MDP(20) + t126 * MDP(21) - t93 * t113 * MDP(22) + (t126 * pkin(10) + t93 * t99) * MDP(23) + (-pkin(9) * MDP(12) - t140 + MDP(9) + (t137 * pkin(10) - MDP(16)) * t116) * t117 + (MDP(8) - pkin(9) * MDP(11) + (-t104 + t106) * MDP(14) + (-pkin(9) * MDP(18) + t124) * t116 + (t116 * MDP(13) + pkin(9) * MDP(19) - t125) * t113) * t114; MDP(10) + t104 * MDP(13) + (t142 * pkin(10) ^ 2 + t99 ^ 2) * MDP(23) + 0.2e1 * t142 * MDP(21) * pkin(10) + 0.2e1 * t125 * t116 + 0.2e1 * (t124 + t138) * t113; t123 * t78 + t131 * t77; t123 * t87 + t131 * t86; t98 * MDP(20) + t92 * MDP(22) + (-t89 * pkin(5) + t88 * qJ(6)) * MDP(23) + ((-0.2e1 * pkin(5) - t150) * MDP(20) + t132) * t117 + (t130 * MDP(21) + t121) * t114 + t122; t140 + t116 * MDP(16) + t129 * MDP(21) + (t131 * t113 + t123 * t116) * pkin(10); 0.2e1 * pkin(5) * MDP(20) + (pkin(5) ^ 2 + (qJ(6) ^ 2)) * MDP(23) - t132; t77 * MDP(23); t86 * MDP(23); t117 * MDP(20) + t89 * MDP(23) + t116 * t141; t152 * t113; t134; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
