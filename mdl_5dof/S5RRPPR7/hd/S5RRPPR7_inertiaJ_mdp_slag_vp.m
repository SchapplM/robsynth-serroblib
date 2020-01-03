% Calculate joint inertia matrix for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPPR7_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:36:23
% EndTime: 2019-12-31 19:36:24
% DurationCPUTime: 0.19s
% Computational Cost: add. (278->84), mult. (510->120), div. (0->0), fcn. (498->6), ass. (0->42)
t71 = sin(pkin(8));
t72 = cos(pkin(8));
t74 = sin(qJ(2));
t76 = cos(qJ(2));
t58 = t71 * t74 - t72 * t76;
t59 = t71 * t76 + t72 * t74;
t68 = -t76 * pkin(2) - pkin(1);
t81 = -t59 * qJ(4) + t68;
t52 = t58 * pkin(3) + t81;
t96 = -0.2e1 * t52;
t95 = 0.2e1 * t76;
t92 = -qJ(3) - pkin(6);
t61 = t92 * t74;
t62 = t92 * t76;
t55 = t71 * t61 - t72 * t62;
t51 = -t58 * pkin(4) + t55;
t94 = t51 * t58;
t73 = sin(qJ(5));
t75 = cos(qJ(5));
t93 = t73 * t75;
t91 = t58 * MDP(14);
t90 = t59 * MDP(21);
t67 = -t72 * pkin(2) - pkin(3);
t89 = t67 * MDP(16);
t88 = t73 * MDP(22);
t87 = t75 * MDP(23);
t86 = MDP(18) * t93;
t53 = -t72 * t61 - t71 * t62;
t85 = t53 ^ 2 + t55 ^ 2;
t83 = t75 * MDP(22) - t73 * MDP(23);
t82 = t87 + t88;
t80 = MDP(13) + t83;
t79 = (MDP(19) * t73 + MDP(20) * t75) * t58;
t78 = t75 * MDP(19) - t73 * MDP(20) + t83 * (-pkin(7) + t67);
t70 = t75 ^ 2;
t69 = t73 ^ 2;
t64 = t71 * pkin(2) + qJ(4);
t50 = t59 * pkin(4) + t53;
t49 = (pkin(3) + pkin(7)) * t58 + t81;
t48 = t75 * t49 + t73 * t50;
t47 = -t73 * t49 + t75 * t50;
t1 = [MDP(1) + pkin(1) * MDP(9) * t95 + (t68 ^ 2 + t85) * MDP(12) + t91 * t96 + (t52 ^ 2 + t85) * MDP(16) + (t69 * MDP(17) + 0.2e1 * t86) * t58 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t74 + MDP(5) * t95) * t74 + (MDP(15) * t96 + 0.2e1 * t79 + t90) * t59 + 0.2e1 * (t47 * t59 - t75 * t94) * MDP(22) + 0.2e1 * (-t48 * t59 + t73 * t94) * MDP(23) + 0.2e1 * (MDP(11) + MDP(13)) * (t53 * t59 - t55 * t58); t74 * MDP(6) + t76 * MDP(7) + t53 * MDP(14) + t55 * MDP(15) + (t53 * t67 + t55 * t64) * MDP(16) + t82 * t51 + (-t76 * MDP(10) - t74 * MDP(9)) * pkin(6) + (t67 * MDP(13) + t78) * t59 + (MDP(17) * t93 + (-t69 + t70) * MDP(18) - t80 * t64) * t58 + ((-t58 * t71 - t59 * t72) * MDP(11) + (-t53 * t72 + t55 * t71) * MDP(12)) * pkin(2); -0.2e1 * t86 + t70 * MDP(17) + MDP(8) + (t71 ^ 2 + t72 ^ 2) * MDP(12) * pkin(2) ^ 2 + (0.2e1 * MDP(14) + t89) * t67 + (MDP(16) * t64 + 0.2e1 * MDP(15) + 0.2e1 * t87 + 0.2e1 * t88) * t64; t68 * MDP(12) - t91 + t52 * MDP(16) + (-MDP(15) - t82) * t59; 0; MDP(12) + MDP(16); t53 * MDP(16) + t80 * t59; MDP(14) + t89; 0; MDP(16); t47 * MDP(22) - t48 * MDP(23) + t79 + t90; t78; -t82; t83; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
