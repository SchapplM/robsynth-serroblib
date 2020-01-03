% Calculate joint inertia matrix for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRRP11_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:33
% EndTime: 2019-12-31 18:54:35
% DurationCPUTime: 0.34s
% Computational Cost: add. (463->108), mult. (857->148), div. (0->0), fcn. (869->6), ass. (0->49)
t111 = cos(qJ(3));
t79 = sin(pkin(8));
t80 = cos(pkin(8));
t82 = sin(qJ(3));
t69 = t111 * t79 + t82 * t80;
t116 = 0.2e1 * t69;
t81 = sin(qJ(4));
t77 = t81 ^ 2;
t83 = cos(qJ(4));
t78 = t83 ^ 2;
t107 = t77 + t78;
t115 = t107 * MDP(23);
t114 = MDP(25) * pkin(7) + MDP(23);
t74 = -pkin(2) * t80 - pkin(1);
t113 = 0.2e1 * t74;
t68 = -t111 * t80 + t79 * t82;
t112 = pkin(4) * t68;
t110 = pkin(1) * MDP(7);
t109 = pkin(6) + qJ(2);
t63 = pkin(3) * t68 - pkin(7) * t69 + t74;
t71 = t109 * t79;
t72 = t109 * t80;
t66 = t111 * t72 - t82 * t71;
t59 = t81 * t63 + t83 * t66;
t105 = qJ(5) * t68;
t104 = t79 * MDP(5);
t103 = t80 * MDP(4);
t102 = MDP(16) * t83;
t101 = t68 * MDP(19);
t100 = MDP(20) + MDP(22);
t99 = -MDP(21) + MDP(24);
t98 = -t83 * t63 + t66 * t81;
t97 = -MDP(25) * pkin(4) - MDP(22);
t96 = -pkin(4) * t83 - qJ(5) * t81;
t95 = pkin(4) * t81 - qJ(5) * t83;
t56 = t105 + t59;
t57 = t98 - t112;
t93 = t56 * t81 - t57 * t83;
t92 = MDP(20) - t97;
t70 = -pkin(3) + t96;
t91 = MDP(20) * pkin(3) - MDP(22) * t70;
t90 = -pkin(3) * MDP(21) - t70 * MDP(24);
t89 = MDP(25) * qJ(5) + t99;
t88 = t83 * MDP(17) - t81 * MDP(18);
t87 = t81 * MDP(17) + t83 * MDP(18);
t86 = -MDP(20) * t98 - t59 * MDP(21);
t65 = t111 * t71 + t82 * t72;
t60 = t95 * t69 + t65;
t1 = [MDP(1) + (t56 ^ 2 + t57 ^ 2 + t60 ^ 2) * MDP(25) + (0.2e1 * t103 - 0.2e1 * t104 + t110) * pkin(1) + (MDP(13) * t113 + t101 + (-MDP(9) + t88) * t116) * t68 + 0.2e1 * (-t57 * MDP(22) + t56 * MDP(24) + t86) * t68 + (-t93 * MDP(23) + (t81 * MDP(20) + t83 * MDP(21)) * t65 + (t81 * MDP(22) - t83 * MDP(24)) * t60) * t116 + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t79 ^ 2 + t80 ^ 2) * qJ(2) + (MDP(14) * t113 + (t78 * MDP(15) - 0.2e1 * t81 * t102 + MDP(8)) * t69) * t69; -t103 + t104 - t110 + t93 * MDP(25) + (MDP(14) - t115) * t69 + (t100 * t83 + t99 * t81 + MDP(13)) * t68; t107 * MDP(25) + MDP(7); -t66 * MDP(14) + (-t83 * MDP(22) - t81 * MDP(24) + MDP(25) * t70) * t60 + (-t83 * MDP(20) + t81 * MDP(21) - MDP(13)) * t65 + (-MDP(11) + (-t100 * t81 + t99 * t83) * pkin(7) + t87) * t68 + (MDP(10) + (-t77 + t78) * MDP(16) + t90 * t83 + (MDP(15) * t83 - t91) * t81) * t69 + t114 * (t56 * t83 + t57 * t81); 0; MDP(12) + t77 * MDP(15) + (t107 * pkin(7) ^ 2 + t70 ^ 2) * MDP(25) + 0.2e1 * pkin(7) * t115 + 0.2e1 * t91 * t83 + 0.2e1 * (t90 + t102) * t81; t101 + (-t98 + 0.2e1 * t112) * MDP(22) + (0.2e1 * t105 + t59) * MDP(24) + (-pkin(4) * t57 + qJ(5) * t56) * MDP(25) + (t96 * MDP(23) + t88) * t69 + t86; t89 * t81 + t92 * t83; -t95 * MDP(23) + (-t92 * t81 + t89 * t83) * pkin(7) + t87; MDP(19) + 0.2e1 * pkin(4) * MDP(22) + 0.2e1 * qJ(5) * MDP(24) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(25); t69 * t83 * MDP(23) - t68 * MDP(22) + MDP(25) * t57; -t83 * MDP(25); t114 * t81; t97; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
