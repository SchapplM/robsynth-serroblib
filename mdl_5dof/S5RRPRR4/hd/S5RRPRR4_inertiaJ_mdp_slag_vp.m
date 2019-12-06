% Calculate joint inertia matrix for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRPRR4_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:32:09
% EndTime: 2019-12-05 18:32:09
% DurationCPUTime: 0.16s
% Computational Cost: add. (209->61), mult. (355->78), div. (0->0), fcn. (317->8), ass. (0->39)
t89 = sin(qJ(2));
t109 = pkin(1) * t89;
t92 = cos(qJ(2));
t80 = t92 * pkin(1) + pkin(2);
t85 = sin(pkin(9));
t86 = cos(pkin(9));
t62 = t86 * t109 + t85 * t80;
t60 = pkin(7) + t62;
t88 = sin(qJ(4));
t56 = (-pkin(8) - t60) * t88;
t91 = cos(qJ(4));
t84 = t91 * pkin(8);
t57 = t91 * t60 + t84;
t87 = sin(qJ(5));
t90 = cos(qJ(5));
t112 = (t90 * t56 - t87 * t57) * MDP(20) + (-t87 * t56 - t90 * t57) * MDP(21);
t77 = t85 * pkin(2) + pkin(7);
t68 = (-pkin(8) - t77) * t88;
t69 = t91 * t77 + t84;
t111 = (t90 * t68 - t87 * t69) * MDP(20) + (-t87 * t68 - t90 * t69) * MDP(21);
t70 = t87 * t88 - t90 * t91;
t71 = t87 * t91 + t90 * t88;
t110 = t70 * MDP(20) + t71 * MDP(21);
t99 = -t91 * MDP(13) + t88 * MDP(14);
t108 = t91 * pkin(4);
t104 = t71 * MDP(17) - t70 * MDP(18);
t78 = -t86 * pkin(2) - pkin(3);
t61 = -t85 * t109 + t86 * t80;
t101 = t88 * MDP(10) + t91 * MDP(11) + t104;
t59 = -pkin(3) - t61;
t100 = MDP(4) + (MDP(8) * t88 + 0.2e1 * MDP(9) * t91) * t88 + (MDP(15) * t71 - 0.2e1 * MDP(16) * t70) * t71;
t98 = -MDP(13) * t88 - MDP(14) * t91;
t97 = (t92 * MDP(5) - t89 * MDP(6)) * pkin(1);
t96 = (MDP(20) * t90 - MDP(21) * t87) * pkin(4);
t95 = 0.2e1 * t99;
t94 = -0.2e1 * t110;
t72 = t78 - t108;
t58 = t59 - t108;
t1 = [MDP(1) + (t61 ^ 2 + t62 ^ 2) * MDP(7) + t59 * t95 - t58 * t94 + 0.2e1 * t97 + t100; (t61 * t86 + t62 * t85) * MDP(7) * pkin(2) + t97 + t100 + t110 * (t58 + t72) + t99 * (t59 + t78); (t85 ^ 2 + t86 ^ 2) * MDP(7) * pkin(2) ^ 2 + t78 * t95 - t72 * t94 + t100; 0; 0; MDP(7); t98 * t60 + t101 + t112; t98 * t77 + t101 + t111; -t99 - t110; MDP(12) + MDP(19) + 0.2e1 * t96; t104 + t112; t104 + t111; -t110; MDP(19) + t96; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
