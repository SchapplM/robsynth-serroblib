% Calculate joint inertia matrix for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5PRRRR6_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:13
% EndTime: 2019-12-05 17:10:14
% DurationCPUTime: 0.21s
% Computational Cost: add. (179->57), mult. (336->73), div. (0->0), fcn. (335->8), ass. (0->35)
t83 = sin(qJ(4));
t87 = cos(qJ(4));
t95 = t87 * MDP(13) - t83 * MDP(14);
t82 = sin(qJ(5));
t86 = cos(qJ(5));
t66 = t82 * t83 - t86 * t87;
t68 = t82 * t87 + t86 * t83;
t107 = t66 * MDP(20) + t68 * MDP(21);
t84 = sin(qJ(3));
t75 = t84 * pkin(2) + pkin(7);
t64 = (-pkin(8) - t75) * t83;
t81 = t87 * pkin(8);
t65 = t87 * t75 + t81;
t109 = (t86 * t64 - t82 * t65) * MDP(20) + (-t82 * t64 - t86 * t65) * MDP(21);
t71 = (-pkin(7) - pkin(8)) * t83;
t72 = t87 * pkin(7) + t81;
t108 = (t86 * t71 - t82 * t72) * MDP(20) + (-t82 * t71 - t86 * t72) * MDP(21);
t88 = cos(qJ(3));
t106 = t88 * pkin(2);
t85 = sin(qJ(2));
t89 = cos(qJ(2));
t69 = t84 * t89 + t88 * t85;
t104 = (-MDP(20) * t68 + MDP(21) * t66) * t69;
t103 = t68 * MDP(17) - t66 * MDP(18);
t77 = -t87 * pkin(4) - pkin(3);
t97 = t83 * MDP(10) + t87 * MDP(11) + t103;
t96 = MDP(5) + (MDP(8) * t83 + 0.2e1 * MDP(9) * t87) * t83 + (MDP(15) * t68 - 0.2e1 * MDP(16) * t66) * t68;
t94 = -MDP(13) * t83 - MDP(14) * t87;
t93 = (t88 * MDP(6) - t84 * MDP(7)) * pkin(2);
t92 = (MDP(20) * t86 - MDP(21) * t82) * pkin(4);
t91 = 0.2e1 * t107;
t90 = -t69 * MDP(7) + (-MDP(6) + t107 - t95) * (t84 * t85 - t88 * t89);
t76 = -pkin(3) - t106;
t70 = t77 - t106;
t1 = [MDP(1); t89 * MDP(3) - t85 * MDP(4) + t90; t70 * t91 - 0.2e1 * t76 * t95 + MDP(2) + 0.2e1 * t93 + t96; t90; t93 + t96 + t95 * (pkin(3) - t76) + t107 * (t70 + t77); 0.2e1 * pkin(3) * t95 + t77 * t91 + t96; t94 * t69 + t104; t94 * t75 + t109 + t97; t94 * pkin(7) + t108 + t97; MDP(12) + MDP(19) + 0.2e1 * t92; t104; t103 + t109; t103 + t108; MDP(19) + t92; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
