% Calculate joint inertia matrix for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(1,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_inertiaJ_mdp_slag_vp: pkin has to be [1x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RPRRR1_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:09:53
% EndTime: 2019-12-05 18:09:54
% DurationCPUTime: 0.24s
% Computational Cost: add. (114->75), mult. (335->116), div. (0->0), fcn. (294->6), ass. (0->37)
t63 = cos(qJ(3));
t58 = sin(qJ(5));
t61 = cos(qJ(5));
t82 = t63 * t61;
t60 = sin(qJ(3));
t62 = cos(qJ(4));
t84 = t60 * t62;
t50 = t58 * t84 + t82;
t83 = t63 * t58;
t51 = t61 * t84 - t83;
t85 = t60 * t61;
t87 = t58 * t60;
t64 = -((-t62 * t83 + t85) * MDP(26) - (t62 * t82 + t87) * MDP(27)) * qJ(2) - t51 * MDP(23) + t50 * MDP(24);
t88 = -t63 * MDP(17) + t64;
t75 = 0.2e1 * t60;
t59 = sin(qJ(4));
t86 = t59 * t62;
t54 = t60 ^ 2;
t57 = t63 ^ 2;
t81 = t54 + t57;
t80 = qJ(2) * t63;
t79 = t51 * MDP(21);
t78 = t59 * MDP(20);
t77 = t61 * MDP(21);
t74 = t58 * t61 * MDP(22);
t72 = t61 * MDP(23) - t58 * MDP(24);
t71 = t58 * MDP(23) + t61 * MDP(24);
t69 = MDP(26) * t58 + MDP(27) * t61;
t68 = -t62 * MDP(19) - MDP(12) + t78;
t67 = -MDP(15) + t72;
t66 = -MDP(17) + t71;
t65 = MDP(26) * t61 - MDP(27) * t58 + MDP(19);
t56 = t62 ^ 2;
t55 = t61 ^ 2;
t53 = t59 ^ 2;
t52 = t58 ^ 2;
t1 = [t57 * MDP(18) + MDP(1) + (-t62 * MDP(16) + MDP(8)) * t63 * t75 + (-0.2e1 * t50 * MDP(22) + t79) * t51 + (t56 * MDP(14) + t53 * MDP(25) + MDP(7)) * t54 + (0.2e1 * t81 * MDP(20) * t62 + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + (-0.2e1 * t62 * t54 * MDP(15) + 0.2e1 * (t81 * MDP(19) + (t50 * MDP(26) + t51 * MDP(27)) * t63) * qJ(2) - t88 * t75) * t59; -MDP(4) + t60 * MDP(13) + (-t62 * t50 - t53 * t87) * MDP(26) + (-t51 * t62 - t53 * t85) * MDP(27) + t68 * t63; MDP(6); (MDP(10) + (t69 * t53 - MDP(13)) * qJ(2)) * t63 + t88 * t62 + (-t63 * MDP(16) + t51 * t77 + (-t50 * t61 - t51 * t58) * MDP(22)) * t59 + (t56 * MDP(15) + MDP(9) + (MDP(14) - MDP(25)) * t86 + t67 * t53 + t68 * qJ(2)) * t60; 0; t56 * MDP(25) + MDP(11) - 0.2e1 * t67 * t86 + (MDP(21) * t55 + MDP(14) - 0.2e1 * t74) * t53; -t63 * MDP(18) + t58 * t79 + (-t58 * t50 + t51 * t61) * MDP(22) + (t60 * MDP(16) - MDP(20) * t80) * t62 + (t66 * t60 - t65 * t80) * t59; t65 * t62 - t78; -t66 * t62 + (MDP(16) + t58 * t77 + (-t52 + t55) * MDP(22)) * t59; t52 * MDP(21) + MDP(18) + 0.2e1 * t74; t59 * t60 * MDP(25) - t64; -t69 * t59; -t62 * MDP(25) + t72 * t59; t71; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
