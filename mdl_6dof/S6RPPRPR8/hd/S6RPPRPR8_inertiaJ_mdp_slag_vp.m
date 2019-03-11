% Calculate joint inertia matrix for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPPRPR8_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:56:12
% EndTime: 2019-03-09 01:56:12
% DurationCPUTime: 0.28s
% Computational Cost: add. (361->95), mult. (588->123), div. (0->0), fcn. (592->6), ass. (0->52)
t89 = sin(qJ(6));
t91 = cos(qJ(6));
t97 = MDP(27) * t91 - t89 * MDP(28);
t95 = -MDP(18) - t97;
t117 = cos(qJ(4));
t86 = sin(pkin(9));
t87 = cos(pkin(9));
t90 = sin(qJ(4));
t71 = -t117 * t87 + t86 * t90;
t67 = t71 ^ 2;
t72 = t117 * t86 + t90 * t87;
t68 = t72 ^ 2;
t122 = t67 + t68;
t119 = pkin(4) + pkin(8);
t121 = -(-MDP(28) * t119 + MDP(25)) * t89 + (-MDP(27) * t119 + MDP(24)) * t91;
t120 = 2 * MDP(20);
t88 = -pkin(1) - qJ(3);
t118 = -pkin(7) + t88;
t75 = t118 * t86;
t76 = t118 * t87;
t65 = t117 * t75 + t90 * t76;
t62 = -t72 * pkin(5) + t65;
t116 = t62 * t72;
t115 = t89 * t91;
t77 = t86 ^ 2 + t87 ^ 2;
t114 = MDP(21) * pkin(4);
t78 = t86 * pkin(3) + qJ(2);
t113 = MDP(27) * t89;
t111 = MDP(28) * t91;
t109 = (MDP(21) * qJ(5));
t108 = MDP(20) - MDP(17);
t107 = 0.2e1 * t72;
t106 = MDP(23) * t115;
t105 = MDP(19) - t114;
t102 = qJ(5) * t71 + t78;
t64 = -t117 * t76 + t90 * t75;
t101 = t64 * t71 + t65 * t72;
t100 = -MDP(16) + t105;
t99 = t86 * MDP(7) + t87 * MDP(8);
t98 = MDP(24) * t89 + MDP(25) * t91;
t96 = t111 + t113;
t94 = t96 + t108;
t93 = qJ(2) ^ 2;
t85 = t91 ^ 2;
t84 = t89 ^ 2;
t70 = t77 * t88;
t63 = pkin(4) * t72 + t102;
t61 = -t71 * pkin(5) + t64;
t60 = t119 * t72 + t102;
t59 = t60 * t91 + t61 * t89;
t58 = -t60 * t89 + t61 * t91;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2) + t93) * MDP(6) + (t77 * t88 ^ 2 + t93) * MDP(10) + (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) * MDP(21) + (MDP(16) * t78 - MDP(19) * t63) * t107 + (MDP(11) + MDP(26)) * t67 + (MDP(22) * t84 + 0.2e1 * t106) * t68 + 0.2e1 * (MDP(5) + t99) * qJ(2) + (-0.2e1 * t78 * MDP(17) + t63 * t120 + (MDP(12) - t98) * t107) * t71 - 0.2e1 * t70 * MDP(9) - 0.2e1 * t101 * MDP(18) + 0.2e1 * (-t91 * t116 - t58 * t71) * MDP(27) + 0.2e1 * (t89 * t116 + t59 * t71) * MDP(28); t70 * MDP(10) + t101 * MDP(21) - pkin(1) * MDP(6) - t77 * MDP(9) + t95 * t122 + MDP(4); t77 * MDP(10) + t122 * MDP(21) + MDP(6); qJ(2) * MDP(10) + MDP(21) * t63 + (MDP(16) - MDP(19)) * t72 + t94 * t71 + t99; 0; MDP(10) + MDP(21); (t108 + t109) * t65 + t100 * t64 + t96 * t62 + (pkin(4) * MDP(18) - MDP(13) - t121) * t71 + (-MDP(14) + MDP(22) * t115 + (-t84 + t85) * MDP(23) + t95 * qJ(5)) * t72; t100 * t71 + (t94 + t109) * t72; 0; -0.2e1 * t106 + t85 * MDP(22) + MDP(15) + (-0.2e1 * MDP(19) + t114) * pkin(4) + (t120 + t109 + 0.2e1 * t111 + 0.2e1 * t113) * qJ(5); t64 * MDP(21) + t95 * t71; t71 * MDP(21); 0; t105; MDP(21); -MDP(26) * t71 + t58 * MDP(27) - t59 * MDP(28) + t98 * t72; t97 * t71; -t96; t121; t97; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
