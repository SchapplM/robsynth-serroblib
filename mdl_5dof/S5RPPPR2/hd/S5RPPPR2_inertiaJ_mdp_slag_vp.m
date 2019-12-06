% Calculate joint inertia matrix for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPPR2_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:31:33
% EndTime: 2019-12-05 17:31:35
% DurationCPUTime: 0.28s
% Computational Cost: add. (338->104), mult. (713->161), div. (0->0), fcn. (712->8), ass. (0->49)
t119 = pkin(1) * MDP(7);
t96 = sin(pkin(7));
t98 = cos(pkin(8));
t116 = t96 * t98;
t94 = sin(pkin(9));
t97 = cos(pkin(9));
t99 = cos(pkin(7));
t76 = t116 * t94 + t99 * t97;
t118 = t94 * t76;
t95 = sin(pkin(8));
t117 = t95 * t96;
t114 = qJ(2) * t99;
t81 = -t99 * pkin(2) - t96 * qJ(3) - pkin(1);
t74 = t98 * t114 + t95 * t81;
t70 = -t99 * qJ(4) + t74;
t75 = (pkin(3) * t95 - qJ(4) * t98 + qJ(2)) * t96;
t66 = t97 * t70 + t94 * t75;
t115 = t94 ^ 2 + t97 ^ 2;
t100 = sin(qJ(5));
t113 = t100 * t95;
t101 = cos(qJ(5));
t112 = t101 * t95;
t111 = t98 * MDP(9);
t110 = t99 * MDP(4);
t77 = t116 * t97 - t99 * t94;
t68 = t77 * t100 - t112 * t96;
t109 = t68 * MDP(19);
t108 = t76 * MDP(12);
t107 = t76 * MDP(20);
t106 = t77 * MDP(13);
t73 = -t95 * t114 + t98 * t81;
t71 = t99 * pkin(3) - t73;
t65 = -t94 * t70 + t97 * t75;
t105 = t71 * MDP(15) + t106;
t104 = (-t101 * t98 - t113 * t97) * MDP(21) - (-t100 * t98 + t112 * t97) * MDP(22);
t103 = t101 * MDP(21) - t100 * MDP(22);
t102 = qJ(2) ^ 2;
t93 = t99 ^ 2;
t92 = t98 ^ 2;
t90 = t96 ^ 2;
t89 = t95 ^ 2;
t86 = t90 * t102;
t69 = t77 * t101 + t113 * t96;
t64 = pkin(6) * t117 + t66;
t63 = -pkin(4) * t117 - t65;
t62 = t76 * pkin(4) - t77 * pkin(6) + t71;
t61 = t100 * t62 + t101 * t64;
t60 = -t100 * t64 + t101 * t62;
t1 = [MDP(1) + (t93 * t102 + t86) * MDP(7) + (t73 ^ 2 + t74 ^ 2 + t86) * MDP(11) + (t65 ^ 2 + t66 ^ 2 + t71 ^ 2) * MDP(15) + (t107 - 0.2e1 * t109) * t76 + (-0.2e1 * t96 * MDP(5) + 0.2e1 * t110 + t119) * pkin(1) + (MDP(16) * t69 - 0.2e1 * t68 * MDP(17) + 0.2e1 * t76 * MDP(18)) * t69 + 0.2e1 * (-t65 * t77 - t66 * t76) * MDP(14) + 0.2e1 * (t60 * t76 + t63 * t68) * MDP(21) + 0.2e1 * (-t61 * t76 + t63 * t69) * MDP(22) + 0.2e1 * (-t73 * MDP(8) + t74 * MDP(9)) * t99 + 0.2e1 * (t106 + t108) * t71 + 0.2e1 * (-t73 * t98 * MDP(10) + (-t74 * MDP(10) + t65 * MDP(12) - t66 * MDP(13)) * t95) * t96 + 0.2e1 * (t93 * MDP(6) + (t95 * MDP(8) + MDP(6) + t111) * t90) * qJ(2); -t110 - t119 + t104 * t76 + (t73 * MDP(11) - t99 * MDP(8) - t105 - t108) * t98 + (-t92 * MDP(10) + MDP(5) + (-t94 * MDP(12) - t97 * MDP(13) - MDP(10)) * t89) * t96 + (t74 * MDP(11) + t99 * MDP(9) + (-t76 * MDP(14) + t66 * MDP(15)) * t97 + (t77 * MDP(14) - t65 * MDP(15) + t68 * MDP(21) + t69 * MDP(22)) * t94) * t95; MDP(7) + (t89 + t92) * MDP(11) + (t115 * t89 + t92) * MDP(15); (-t97 * t77 - t118) * MDP(14) + (t65 * t97 + t66 * t94) * MDP(15) + (-t100 * t118 - t97 * t68) * MDP(21) + (-t101 * t118 - t97 * t69) * MDP(22) + (qJ(2) * MDP(11) + t111 + (MDP(12) * t97 - MDP(13) * t94 + MDP(8)) * t95) * t96; 0; MDP(15) * t115 + MDP(11); (MDP(12) + t103) * t76 + t105; -t98 * MDP(15); 0; MDP(15); t69 * MDP(18) + t60 * MDP(21) - t61 * MDP(22) + t107 - t109; t104; (-MDP(21) * t100 - MDP(22) * t101) * t94; t103; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
