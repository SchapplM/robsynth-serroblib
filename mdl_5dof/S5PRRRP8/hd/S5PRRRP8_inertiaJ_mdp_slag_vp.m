% Calculate joint inertia matrix for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP8_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:37
% EndTime: 2019-12-05 17:00:38
% DurationCPUTime: 0.33s
% Computational Cost: add. (282->116), mult. (606->172), div. (0->0), fcn. (568->8), ass. (0->48)
t79 = sin(qJ(3));
t112 = 0.2e1 * t79;
t111 = MDP(22) * pkin(8) + MDP(20);
t110 = MDP(17) + MDP(19);
t78 = sin(qJ(4));
t109 = pkin(7) * t78;
t81 = cos(qJ(4));
t108 = pkin(7) * t81;
t76 = sin(pkin(5));
t80 = sin(qJ(2));
t106 = t76 * t80;
t77 = cos(pkin(5));
t82 = cos(qJ(3));
t65 = t79 * t106 - t77 * t82;
t107 = t65 * t79;
t83 = cos(qJ(2));
t105 = t76 * t83;
t104 = t78 * t82;
t70 = -t82 * pkin(3) - t79 * pkin(8) - pkin(2);
t63 = t82 * t108 + t78 * t70;
t73 = t78 ^ 2;
t75 = t81 ^ 2;
t103 = t73 + t75;
t101 = MDP(11) * t79;
t100 = MDP(20) * t79;
t99 = t78 * MDP(14);
t98 = t81 * MDP(13);
t97 = (MDP(18) - MDP(21));
t95 = -MDP(22) * pkin(4) - MDP(19);
t94 = -2 * qJ(5) * MDP(21) - MDP(16);
t93 = -t81 * pkin(4) - t78 * qJ(5);
t92 = -pkin(4) * t78 + t81 * qJ(5);
t59 = -t82 * qJ(5) + t63;
t68 = t81 * t70;
t60 = -t68 + (pkin(4) + t109) * t82;
t90 = t59 * t81 + t60 * t78;
t89 = -MDP(17) + t95;
t88 = MDP(22) * qJ(5) - t97;
t87 = t81 * MDP(14) - t78 * MDP(15);
t86 = (-pkin(7) * t104 + t68) * MDP(17) - t63 * MDP(18);
t85 = t78 * MDP(19) - t81 * MDP(21);
t71 = pkin(8) * t104;
t69 = -pkin(3) + t93;
t66 = t82 * t106 + t77 * t79;
t64 = (pkin(7) - t92) * t79;
t58 = -t78 * t105 + t66 * t81;
t57 = t81 * t105 + t66 * t78;
t1 = [MDP(1) + (t57 ^ 2 + t58 ^ 2 + t65 ^ 2) * MDP(22); (t57 * t60 + t58 * t59 + t65 * t64) * MDP(22) + (t57 * t81 - t58 * t78) * t100 + t97 * (t81 * t107 + t58 * t82) + (-t80 * MDP(4) + (MDP(10) * t82 + MDP(3) - t101) * t83) * t76 + t110 * (t78 * t107 + t57 * t82); MDP(2) - 0.2e1 * pkin(2) * t101 + (t59 ^ 2 + t60 ^ 2 + t64 ^ 2) * MDP(22) + (0.2e1 * pkin(2) * MDP(10) + t82 * MDP(16) + (MDP(6) - t87) * t112) * t82 + 0.2e1 * (t60 * MDP(19) - t59 * MDP(21) - t86) * t82 + ((-t59 * t78 + t60 * t81) * MDP(20) + t85 * t64) * t112 + (t75 * MDP(12) - 0.2e1 * t78 * t98 + MDP(5) + 0.2e1 * (t78 * MDP(17) + t81 * MDP(18)) * pkin(7)) * t79 ^ 2; -t66 * MDP(11) + (t69 * MDP(22) - t110 * t81 + t97 * t78 - MDP(10)) * t65 + t111 * (t57 * t78 + t58 * t81); t71 * MDP(17) + (-t64 * t81 + t71) * MDP(19) + t90 * MDP(20) - t64 * t78 * MDP(21) + (t90 * pkin(8) + t64 * t69) * MDP(22) + (-pkin(7) * MDP(11) - t99 + MDP(8) + (t97 * pkin(8) - MDP(15)) * t81) * t82 + (MDP(7) - pkin(7) * MDP(10) + t81 * t78 * MDP(12) + (-t73 + t75) * MDP(13) + (-pkin(3) * t78 - t108) * MDP(17) + (-pkin(3) * t81 + t109) * MDP(18) + t85 * t69) * t79; MDP(9) + t73 * MDP(12) + (t103 * pkin(8) ^ 2 + t69 ^ 2) * MDP(22) + 0.2e1 * t103 * MDP(20) * pkin(8) + 0.2e1 * (pkin(3) * MDP(17) - t69 * MDP(19)) * t81 + 0.2e1 * (-pkin(3) * MDP(18) - t69 * MDP(21) + t98) * t78; t89 * t57 + t88 * t58; t68 * MDP(19) + t63 * MDP(21) + (-t60 * pkin(4) + t59 * qJ(5)) * MDP(22) + ((-0.2e1 * pkin(4) - t109) * MDP(19) + t94) * t82 + (t93 * MDP(20) + t87) * t79 + t86; t99 + t81 * MDP(15) + t92 * MDP(20) + (t89 * t78 + t88 * t81) * pkin(8); 0.2e1 * pkin(4) * MDP(19) + (pkin(4) ^ 2 + (qJ(5) ^ 2)) * MDP(22) - t94; t57 * MDP(22); t82 * MDP(19) + t60 * MDP(22) + t81 * t100; t111 * t78; t95; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
