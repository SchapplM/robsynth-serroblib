% Calculate joint inertia matrix for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPRRR5_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:29
% EndTime: 2019-12-05 18:16:30
% DurationCPUTime: 0.18s
% Computational Cost: add. (187->56), mult. (327->70), div. (0->0), fcn. (293->8), ass. (0->37)
t76 = cos(pkin(9));
t68 = pkin(1) * t76 + pkin(2);
t79 = sin(qJ(3));
t82 = cos(qJ(3));
t75 = sin(pkin(9));
t99 = pkin(1) * t75;
t56 = -t68 * t79 - t82 * t99;
t54 = pkin(7) - t56;
t78 = sin(qJ(4));
t50 = (-pkin(8) - t54) * t78;
t81 = cos(qJ(4));
t74 = t81 * pkin(8);
t51 = t54 * t81 + t74;
t77 = sin(qJ(5));
t80 = cos(qJ(5));
t102 = (t50 * t80 - t51 * t77) * MDP(20) + (-t50 * t77 - t51 * t80) * MDP(21);
t64 = (-pkin(7) - pkin(8)) * t78;
t65 = pkin(7) * t81 + t74;
t101 = (t64 * t80 - t65 * t77) * MDP(20) + (-t64 * t77 - t65 * t80) * MDP(21);
t62 = t77 * t78 - t80 * t81;
t63 = t77 * t81 + t78 * t80;
t100 = t62 * MDP(20) + t63 * MDP(21);
t87 = MDP(13) * t81 - MDP(14) * t78;
t98 = pkin(4) * t81;
t94 = t63 * MDP(17) - t62 * MDP(18);
t55 = t68 * t82 - t79 * t99;
t93 = t55 * MDP(6);
t92 = t56 * MDP(7);
t89 = t78 * MDP(10) + t81 * MDP(11) + t94;
t53 = -pkin(3) - t55;
t88 = MDP(5) + (MDP(8) * t78 + 0.2e1 * MDP(9) * t81) * t78 + (MDP(15) * t63 - 0.2e1 * MDP(16) * t62) * t63;
t86 = -MDP(13) * t78 - MDP(14) * t81;
t85 = (MDP(20) * t80 - MDP(21) * t77) * pkin(4);
t84 = -0.2e1 * t100;
t70 = -pkin(3) - t98;
t52 = t53 - t98;
t1 = [MDP(1) + (t75 ^ 2 + t76 ^ 2) * MDP(4) * pkin(1) ^ 2 - 0.2e1 * t87 * t53 - t52 * t84 + 0.2e1 * t93 + 0.2e1 * t92 + t88; 0; MDP(4); t88 + t92 + t93 + t87 * (pkin(3) - t53) + t100 * (t52 + t70); 0; 0.2e1 * pkin(3) * t87 - t70 * t84 + t88; t86 * t54 + t102 + t89; t87 - t100; pkin(7) * t86 + t101 + t89; MDP(12) + MDP(19) + 0.2e1 * t85; t94 + t102; -t100; t94 + t101; MDP(19) + t85; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
