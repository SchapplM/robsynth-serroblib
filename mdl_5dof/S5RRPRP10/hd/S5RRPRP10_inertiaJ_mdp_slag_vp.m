% Calculate joint inertia matrix for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRP10_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:07
% EndTime: 2019-12-31 20:11:08
% DurationCPUTime: 0.28s
% Computational Cost: add. (228->94), mult. (388->119), div. (0->0), fcn. (290->4), ass. (0->45)
t97 = pkin(3) + pkin(6);
t73 = -pkin(2) - pkin(7);
t70 = sin(qJ(2));
t59 = t97 * t70;
t69 = sin(qJ(4));
t95 = t69 * t59;
t71 = cos(qJ(4));
t94 = t69 * t71;
t72 = cos(qJ(2));
t60 = t97 * t72;
t93 = MDP(22) * pkin(4);
t92 = MDP(23) * pkin(4);
t91 = pkin(2) * MDP(14);
t90 = -qJ(5) + t73;
t89 = MDP(18) * t71;
t88 = t69 * MDP(20);
t87 = t69 * MDP(21);
t86 = t71 * MDP(21);
t85 = pkin(6) ^ 2 * MDP(14);
t84 = MDP(16) * t94;
t83 = -t70 * qJ(3) - pkin(1);
t54 = t73 * t72 + t83;
t82 = qJ(5) * t72 - t54;
t81 = MDP(12) - t91;
t80 = MDP(20) * t73 + MDP(17);
t55 = t71 * t59;
t48 = t70 * pkin(4) + t82 * t69 + t55;
t49 = -t82 * t71 + t95;
t79 = t48 * t71 + t49 * t69;
t78 = (-MDP(21) * t73 - MDP(18)) * t69;
t77 = (-t69 * t54 + t55) * MDP(20) - (t71 * t54 + t95) * MDP(21);
t76 = t71 * MDP(20) - t87;
t75 = pkin(6) * MDP(14) + MDP(11) + t76;
t68 = t72 ^ 2;
t67 = t71 ^ 2;
t66 = t70 ^ 2;
t65 = t69 ^ 2;
t62 = t69 * pkin(4) + qJ(3);
t61 = t65 + t67;
t58 = -t72 * pkin(2) + t83;
t57 = t90 * t71;
t56 = t90 * t69;
t53 = t71 * t72 * pkin(4) + t60;
t52 = t56 * t69 + t57 * t71;
t1 = [MDP(1) + (t48 ^ 2 + t49 ^ 2 + t53 ^ 2) * MDP(23) + t58 ^ 2 * MDP(14) + (t65 * MDP(15) + 0.2e1 * t84 + t85) * t68 + (MDP(19) + MDP(4) + t85) * t66 + 0.2e1 * (t66 + t68) * MDP(11) * pkin(6) + 0.2e1 * (pkin(1) * MDP(9) + t58 * MDP(12) + (t48 * t69 - t49 * t71) * MDP(22) + t76 * t60) * t72 + 0.2e1 * (-pkin(1) * MDP(10) - t58 * MDP(13) + (-MDP(17) * t69 + MDP(5) - t89) * t72 + t77) * t70; -t79 * MDP(22) + (t48 * t57 + t49 * t56 + t53 * t62) * MDP(23) + (t86 + t88) * t60 + (-pkin(2) * MDP(11) + MDP(6) + t80 * t71 + t78 + (-MDP(9) + t81) * pkin(6)) * t70 + (MDP(7) - MDP(15) * t94 + (t65 - t67) * MDP(16) + (-t56 * t71 + t57 * t69) * MDP(22) + (-MDP(10) + MDP(13)) * pkin(6) + t75 * qJ(3)) * t72; MDP(8) + t67 * MDP(15) - 0.2e1 * t84 - 0.2e1 * t52 * MDP(22) + (t56 ^ 2 + t57 ^ 2 + t62 ^ 2) * MDP(23) + (-0.2e1 * MDP(12) + t91) * pkin(2) + (MDP(14) * qJ(3) + 0.2e1 * MDP(13) + 0.2e1 * t86 + 0.2e1 * t88) * qJ(3); t79 * MDP(23) + t75 * t70; -t61 * MDP(22) + t52 * MDP(23) + t81; t61 * MDP(23) + MDP(14); t48 * t92 + t70 * MDP(19) + (-t89 + (-MDP(17) + t93) * t69) * t72 + t77; t57 * t92 + t78 + (t80 - t93) * t71; -t87 + (MDP(20) + t92) * t71; MDP(23) * pkin(4) ^ 2 + MDP(19); t53 * MDP(23); t62 * MDP(23); 0; 0; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
