% Calculate joint inertia matrix for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP4_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:56:37
% EndTime: 2021-01-15 12:56:40
% DurationCPUTime: 0.43s
% Computational Cost: add. (425->102), mult. (843->136), div. (0->0), fcn. (806->6), ass. (0->46)
t84 = sin(qJ(3));
t86 = cos(qJ(3));
t83 = sin(qJ(4));
t85 = cos(qJ(4));
t70 = -t83 * t84 + t85 * t86;
t71 = t83 * t86 + t84 * t85;
t97 = MDP(20) + MDP(22);
t95 = -t97 * t71 + (MDP(19) + MDP(21)) * t70;
t117 = -t86 * MDP(12) + t84 * MDP(13) - t95;
t103 = qJ(2) * t84;
t81 = sin(pkin(8));
t82 = cos(pkin(8));
t72 = -pkin(2) * t82 - pkin(6) * t81 - pkin(1);
t68 = t86 * t72;
t96 = t86 * t82 * qJ(2);
t116 = (-t82 * t103 + t68) * MDP(12) - (t72 * t84 + t96) * MDP(13) - (MDP(10) * t84 - MDP(9) * t86) * t81;
t62 = t71 * t81;
t63 = t70 * t81;
t111 = t63 * MDP(16) - t62 * MDP(17);
t110 = 0.2e1 * MDP(21);
t109 = pkin(3) * t83;
t108 = pkin(7) * t81;
t105 = (pkin(3) * t84 + qJ(2)) * t81;
t104 = MDP(24) * pkin(4);
t102 = qJ(2) ^ 2 * MDP(6);
t78 = t85 * pkin(3);
t76 = t78 + pkin(4);
t101 = MDP(24) * t76;
t55 = -t86 * t108 + t68 + (-pkin(3) - t103) * t82;
t57 = t96 + (t72 - t108) * t84;
t52 = t85 * t55 - t83 * t57;
t100 = t52 * MDP(19);
t99 = t85 * MDP(19);
t98 = MDP(11) + MDP(18);
t53 = t55 * t83 + t57 * t85;
t91 = MDP(21) * t62 + MDP(22) * t63;
t90 = -t63 * qJ(5) + t52;
t51 = -qJ(5) * t62 + t53;
t89 = -t53 * MDP(20) - t51 * MDP(22) + t100;
t88 = t90 * MDP(21) + t111;
t80 = t82 ^ 2;
t79 = t81 ^ 2;
t74 = t82 * t109;
t56 = pkin(4) * t62 + t105;
t50 = -pkin(4) * t82 + t90;
t1 = [MDP(1) + pkin(1) ^ 2 * MDP(6) + (t50 ^ 2 + t51 ^ 2 + t56 ^ 2) * MDP(24) + (t98 + t102) * t80 + (t102 + (MDP(7) * t86 - 0.2e1 * MDP(8) * t84) * t86) * t79 + (MDP(14) * t63 - 0.2e1 * MDP(15) * t62) * t63 + 0.2e1 * (-t50 * t63 - t51 * t62) * MDP(23) + 0.2e1 * (MDP(19) * t62 + MDP(20) * t63) * t105 + 0.2e1 * t91 * t56 + 0.2e1 * (t80 * MDP(5) + (MDP(12) * t84 + MDP(13) * t86 + MDP(5)) * t79) * qJ(2) + 0.2e1 * (-MDP(21) * t50 + pkin(1) * MDP(4) - t111 - t116 - t89) * t82; -pkin(1) * MDP(6) + (-t62 * t71 - t63 * t70) * MDP(23) + (t50 * t70 + t51 * t71) * MDP(24) + (-MDP(4) + t117) * t82; MDP(6) + (t70 ^ 2 + t71 ^ 2) * MDP(24); t100 + (-t53 + t74) * MDP(20) + (-t51 + t74) * MDP(22) + (-t62 * t109 - t63 * t76) * MDP(23) + (t51 * t109 + t50 * t76) * MDP(24) + (-pkin(3) * t99 + (-pkin(4) - t76) * MDP(21) - t98) * t82 + t88 + t116; (t71 * t109 + t70 * t76) * MDP(24) - t117; (t110 + t101) * t76 + (0.2e1 * t99 + (MDP(24) * t109 - 0.2e1 * MDP(20) - 0.2e1 * MDP(22)) * t83) * pkin(3) + t98; -t82 * MDP(18) + (-0.2e1 * MDP(21) * t82 - MDP(23) * t63 + MDP(24) * t50) * pkin(4) + t88 + t89; t70 * t104 + t95; MDP(18) + (0.2e1 * pkin(4) + t78) * MDP(21) + pkin(4) * t101 + (-t97 * t83 + t99) * pkin(3); MDP(18) + (t110 + t104) * pkin(4); MDP(24) * t56 + t91; 0; 0; 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
