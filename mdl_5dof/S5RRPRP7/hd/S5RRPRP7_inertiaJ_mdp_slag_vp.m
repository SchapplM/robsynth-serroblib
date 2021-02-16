% Calculate joint inertia matrix for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP7_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:41:26
% EndTime: 2021-01-15 20:41:30
% DurationCPUTime: 0.40s
% Computational Cost: add. (517->122), mult. (939->177), div. (0->0), fcn. (939->6), ass. (0->49)
t110 = cos(qJ(2));
t83 = sin(pkin(8));
t84 = cos(pkin(8));
t86 = sin(qJ(2));
t73 = t83 * t110 + t84 * t86;
t115 = 0.2e1 * t73;
t114 = -qJ(3) - pkin(6);
t85 = sin(qJ(4));
t81 = t85 ^ 2;
t87 = cos(qJ(4));
t82 = t87 ^ 2;
t107 = t81 + t82;
t113 = t107 * MDP(23);
t80 = -t110 * pkin(2) - pkin(1);
t112 = 0.2e1 * t80;
t72 = -t84 * t110 + t83 * t86;
t111 = pkin(4) * t72;
t78 = pkin(2) * t83 + pkin(7);
t109 = t72 * t78;
t108 = t73 * t87;
t65 = t72 * pkin(3) - t73 * pkin(7) + t80;
t101 = t114 * t86;
t75 = t114 * t110;
t69 = t83 * t101 - t84 * t75;
t61 = t85 * t65 + t87 * t69;
t106 = qJ(5) * t72;
t105 = MDP(16) * t87;
t104 = t72 * MDP(19);
t103 = 0.2e1 * t110;
t102 = -MDP(21) + MDP(24);
t79 = -pkin(2) * t84 - pkin(3);
t100 = -t87 * t65 + t69 * t85;
t67 = -t84 * t101 - t75 * t83;
t99 = -MDP(25) * pkin(4) - MDP(22);
t98 = t78 * MDP(25) + MDP(23);
t97 = -pkin(4) * t87 - qJ(5) * t85;
t96 = pkin(4) * t85 - qJ(5) * t87;
t58 = t106 + t61;
t59 = t100 - t111;
t95 = t58 * t85 - t59 * t87;
t70 = t79 + t97;
t94 = t70 * t73 - t109;
t93 = t73 * t79 - t109;
t92 = MDP(20) - t99;
t91 = MDP(25) * qJ(5) + t102;
t90 = t87 * MDP(17) - t85 * MDP(18);
t89 = -MDP(20) * t100 - t61 * MDP(21);
t62 = t96 * t73 + t67;
t1 = [MDP(1) + pkin(1) * MDP(9) * t103 + (t67 ^ 2 + t69 ^ 2 + t80 ^ 2) * MDP(14) + (t58 ^ 2 + t59 ^ 2 + t62 ^ 2) * MDP(25) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t86 + MDP(5) * t103) * t86 + (MDP(11) * t112 + t90 * t115 + t104) * t72 + 0.2e1 * (-t69 * MDP(13) - t59 * MDP(22) + t58 * MDP(24) + t89) * t72 + (-t95 * MDP(23) + (t85 * MDP(22) - t87 * MDP(24)) * t62 + (t85 * MDP(20) + t87 * MDP(21) + MDP(13)) * t67) * t115 + (MDP(12) * t112 + (t82 * MDP(15) - 0.2e1 * t85 * t105) * t73) * t73; t62 * t70 * MDP(25) - t67 * MDP(11) - t69 * MDP(12) + t86 * MDP(6) + t110 * MDP(7) + (-t81 + t82) * MDP(16) * t73 + (-t110 * MDP(10) - t86 * MDP(9)) * pkin(6) + (t72 * MDP(18) - t67 * MDP(20) + t93 * MDP(21) - t62 * MDP(22) - t94 * MDP(24) + t98 * t58) * t87 + (MDP(15) * t108 + t72 * MDP(17) + t93 * MDP(20) + t67 * MDP(21) + t94 * MDP(22) - t62 * MDP(24) + t98 * t59) * t85 + ((-t72 * t83 - t73 * t84) * MDP(13) + (-t67 * t84 + t69 * t83) * MDP(14)) * pkin(2); MDP(8) + t81 * MDP(15) + (t107 * t78 ^ 2 + t70 ^ 2) * MDP(25) + (t83 ^ 2 + t84 ^ 2) * MDP(14) * pkin(2) ^ 2 + 0.2e1 * t78 * t113 + 0.2e1 * (-t79 * MDP(20) - t70 * MDP(22)) * t87 + 0.2e1 * (MDP(11) * t84 - MDP(12) * t83) * pkin(2) + 0.2e1 * (MDP(21) * t79 - MDP(24) * t70 + t105) * t85; t80 * MDP(14) + t95 * MDP(25) + (MDP(12) - t113) * t73 + (MDP(11) + (MDP(20) + MDP(22)) * t87 + t102 * t85) * t72; 0; t107 * MDP(25) + MDP(14); t104 + (-t100 + 0.2e1 * t111) * MDP(22) + (0.2e1 * t106 + t61) * MDP(24) + (-pkin(4) * t59 + qJ(5) * t58) * MDP(25) + (t97 * MDP(23) + t90) * t73 + t89; t85 * MDP(17) + t87 * MDP(18) - t96 * MDP(23) + (-t92 * t85 + t91 * t87) * t78; t91 * t85 + t92 * t87; MDP(19) + 0.2e1 * pkin(4) * MDP(22) + 0.2e1 * qJ(5) * MDP(24) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(25); -t72 * MDP(22) + MDP(23) * t108 + t59 * MDP(25); t98 * t85; -t87 * MDP(25); t99; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
