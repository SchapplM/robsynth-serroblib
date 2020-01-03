% Calculate joint inertia matrix for
% S5RRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPRP2_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:43
% EndTime: 2019-12-31 19:49:44
% DurationCPUTime: 0.17s
% Computational Cost: add. (215->65), mult. (349->91), div. (0->0), fcn. (245->6), ass. (0->38)
t73 = sin(pkin(8));
t64 = t73 * pkin(2) + pkin(7);
t75 = sin(qJ(4));
t71 = t75 ^ 2;
t77 = cos(qJ(4));
t92 = t77 ^ 2 + t71;
t93 = t92 * t64;
t99 = 2 * MDP(16);
t76 = sin(qJ(2));
t98 = pkin(1) * t76;
t74 = cos(pkin(8));
t97 = t74 * pkin(2);
t78 = cos(qJ(2));
t66 = t78 * pkin(1) + pkin(2);
t53 = t74 * t66 - t73 * t98;
t83 = -t77 * pkin(4) - t75 * qJ(5) - pkin(3);
t47 = -t53 + t83;
t55 = t83 - t97;
t96 = -t47 - t55;
t54 = t73 * t66 + t74 * t98;
t52 = pkin(7) + t54;
t95 = t92 * t52;
t51 = -pkin(3) - t53;
t65 = -pkin(3) - t97;
t94 = t51 + t65;
t91 = MDP(18) * t75;
t90 = t75 * MDP(10) + t77 * MDP(11) + (-t75 * pkin(4) + t77 * qJ(5)) * MDP(16);
t89 = 0.2e1 * t75 * t77 * MDP(9) + t71 * MDP(8) + MDP(4);
t88 = t92 * MDP(18);
t87 = -MDP(18) * pkin(4) - MDP(15);
t86 = MDP(13) - t87;
t85 = MDP(18) * qJ(5) - MDP(14) + MDP(17);
t84 = -0.2e1 * t77 * MDP(15) - 0.2e1 * t75 * MDP(17);
t82 = (t78 * MDP(5) - t76 * MDP(6)) * pkin(1);
t81 = -0.2e1 * t77 * MDP(13) + 0.2e1 * t75 * MDP(14);
t80 = -t86 * t75 + t85 * t77;
t67 = t75 * MDP(16);
t1 = [MDP(1) + (t53 ^ 2 + t54 ^ 2) * MDP(7) + t95 * t99 + t51 * t81 + t52 ^ 2 * t88 + (t47 * MDP(18) + t84) * t47 + 0.2e1 * t82 + t89; (t93 + t95) * MDP(16) + (t47 * t55 + t52 * t93) * MDP(18) + (t53 * t74 + t54 * t73) * MDP(7) * pkin(2) + t82 + (-t94 * MDP(13) + t96 * MDP(15)) * t77 + (t94 * MDP(14) + t96 * MDP(17)) * t75 + t89; t93 * t99 + (t73 ^ 2 + t74 ^ 2) * MDP(7) * pkin(2) ^ 2 + t65 * t81 + t64 ^ 2 * t88 + (MDP(18) * t55 + t84) * t55 + t89; 0; 0; MDP(7) + t88; t80 * t52 + t90; t80 * t64 + t90; t85 * t75 + t86 * t77; MDP(12) + 0.2e1 * pkin(4) * MDP(15) + 0.2e1 * qJ(5) * MDP(17) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(18); t52 * t91 + t67; t64 * t91 + t67; -t77 * MDP(18); t87; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
