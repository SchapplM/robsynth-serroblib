% Calculate joint inertia matrix for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPPRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPPRR3_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:59
% EndTime: 2019-03-09 01:34:00
% DurationCPUTime: 0.30s
% Computational Cost: add. (344->89), mult. (560->127), div. (0->0), fcn. (563->8), ass. (0->55)
t101 = cos(qJ(6));
t99 = sin(qJ(6));
t107 = MDP(26) * t101 - MDP(27) * t99;
t100 = sin(qJ(5));
t124 = cos(qJ(5));
t95 = sin(pkin(10));
t97 = cos(pkin(10));
t106 = -t100 * t95 + t124 * t97;
t76 = t100 * t97 + t124 * t95;
t71 = t76 * MDP(20);
t129 = MDP(19) * t106 - t71;
t122 = t95 ^ 2 + t97 ^ 2;
t111 = t122 * MDP(13);
t102 = -pkin(1) - pkin(2);
t96 = sin(pkin(9));
t98 = cos(pkin(9));
t82 = t98 * qJ(2) + t96 * t102;
t78 = -qJ(4) + t82;
t128 = t78 * t111;
t108 = MDP(26) * t99 + MDP(27) * t101;
t127 = t99 * MDP(23) + t101 * MDP(24) - t108 * pkin(8);
t80 = t96 * qJ(2) - t98 * t102;
t79 = pkin(3) + t80;
t69 = t97 * pkin(4) + t79;
t126 = -0.2e1 * t69;
t125 = pkin(7) - t78;
t123 = t76 * t99;
t121 = t101 * t76;
t120 = t101 * t99;
t119 = MDP(10) * t97;
t118 = MDP(11) * t95;
t117 = MDP(13) * t79;
t115 = t106 * MDP(25);
t113 = MDP(22) * t120;
t112 = t122 * MDP(12);
t109 = -MDP(23) * t101 + MDP(24) * t99;
t105 = t107 * t106 + t129;
t104 = -MDP(19) - t107;
t103 = t117 - t118 + t119;
t94 = t101 ^ 2;
t93 = t99 ^ 2;
t92 = t98 ^ 2;
t90 = t96 ^ 2;
t67 = t106 * t96;
t66 = t76 * t96;
t65 = t125 * t97;
t64 = t125 * t95;
t63 = t101 * t67 - t98 * t99;
t62 = -t101 * t98 - t67 * t99;
t61 = pkin(5) * t106 + pkin(8) * t76 + t69;
t60 = t100 * t64 - t124 * t65;
t59 = -t100 * t65 - t124 * t64;
t58 = t101 * t60 + t61 * t99;
t57 = t101 * t61 - t60 * t99;
t1 = [MDP(1) + (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + (t80 ^ 2 + t82 ^ 2) * MDP(9) + t71 * t126 + (t117 - 0.2e1 * t118 + 0.2e1 * t119) * t79 + (MDP(21) * t94 + MDP(14) - 0.2e1 * t113) * t76 ^ 2 - (MDP(19) * t126 - t115 + 0.2e1 * (-MDP(15) - t109) * t76) * t106 + 0.2e1 * t80 * MDP(7) + 0.2e1 * t82 * MDP(8) + 0.2e1 * (t106 * t57 - t59 * t123) * MDP(26) + 0.2e1 * (-t106 * t58 - t59 * t121) * MDP(27) + (-0.2e1 * t112 + t128) * t78; -MDP(4) - pkin(1) * MDP(6) + (t106 * t62 - t66 * t123) * MDP(26) + (-t106 * t63 - t66 * t121) * MDP(27) + (t82 * MDP(9) + MDP(8) - t112 + t128) * t96 + (-MDP(9) * t80 - MDP(7) - t103 - t129) * t98; MDP(6) + (t90 + t92) * MDP(9) + (t122 * t90 + t92) * MDP(13); 0; 0; MDP(9) + t111; t103 + t105; -t98 * MDP(13); 0; MDP(13); -t60 * MDP(20) + t104 * t59 - (MDP(17) - t127) * t106 + (-MDP(16) - MDP(21) * t120 + (t93 - t94) * MDP(22) + t108 * pkin(5)) * t76; -MDP(20) * t67 + t104 * t66; t105; 0; MDP(21) * t93 + 0.2e1 * pkin(5) * t107 + MDP(18) + 0.2e1 * t113; t57 * MDP(26) - t58 * MDP(27) + t109 * t76 + t115; MDP(26) * t62 - MDP(27) * t63; -t108 * t76; t107; t127; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
