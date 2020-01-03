% Calculate joint inertia matrix for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR13_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR13_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRPR13_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:54
% EndTime: 2019-12-31 18:32:54
% DurationCPUTime: 0.22s
% Computational Cost: add. (254->81), mult. (474->106), div. (0->0), fcn. (470->6), ass. (0->44)
t97 = -2 * MDP(16);
t96 = pkin(3) + pkin(7);
t95 = cos(qJ(3));
t94 = pkin(1) * MDP(7);
t67 = sin(pkin(8));
t91 = pkin(6) + qJ(2);
t59 = t91 * t67;
t68 = cos(pkin(8));
t60 = t91 * t68;
t70 = sin(qJ(3));
t53 = -t70 * t59 + t95 * t60;
t57 = t70 * t67 - t95 * t68;
t50 = -t57 * pkin(4) + t53;
t93 = t50 * t57;
t69 = sin(qJ(5));
t71 = cos(qJ(5));
t92 = t69 * t71;
t89 = (MDP(18) * pkin(3));
t88 = t67 * MDP(5);
t87 = t68 * MDP(4);
t86 = t69 * MDP(24);
t85 = t71 * MDP(25);
t84 = MDP(18) * qJ(4);
t83 = MDP(17) - MDP(14);
t58 = t95 * t67 + t70 * t68;
t82 = 0.2e1 * t58;
t81 = MDP(20) * t92;
t62 = -t68 * pkin(2) - pkin(1);
t80 = MDP(16) - t89;
t79 = MDP(21) * t69 + MDP(22) * t71;
t78 = t71 * MDP(24) - t69 * MDP(25);
t77 = t85 + t86;
t76 = -t58 * qJ(4) + t62;
t52 = t95 * t59 + t70 * t60;
t75 = MDP(15) + t78;
t74 = (-MDP(24) * t96 + MDP(21)) * t71 + (MDP(25) * t96 - MDP(22)) * t69;
t66 = t71 ^ 2;
t65 = t69 ^ 2;
t51 = t57 * pkin(3) + t76;
t49 = t58 * pkin(4) + t52;
t48 = t96 * t57 + t76;
t47 = t71 * t48 + t69 * t49;
t46 = -t69 * t48 + t71 * t49;
t1 = [MDP(1) + (t51 ^ 2 + t52 ^ 2 + t53 ^ 2) * MDP(18) + (t62 * MDP(14) - t51 * MDP(17)) * t82 + (MDP(8) + MDP(23)) * t58 ^ 2 + (t65 * MDP(19) + 0.2e1 * t81) * t57 ^ 2 + (0.2e1 * t87 - 0.2e1 * t88 + t94) * pkin(1) + (0.2e1 * t62 * MDP(13) + t51 * t97 + (-MDP(9) + t79) * t82) * t57 + 0.2e1 * (t52 * t58 - t53 * t57) * MDP(15) + 0.2e1 * (t46 * t58 - t71 * t93) * MDP(24) + 0.2e1 * (-t47 * t58 + t69 * t93) * MDP(25) + (MDP(7) * qJ(2) + 2 * MDP(6)) * (t67 ^ 2 + t68 ^ 2) * qJ(2); t51 * MDP(18) - t87 + t88 - t94 + (MDP(13) - MDP(16)) * t57 + (-t77 - t83) * t58; MDP(7) + MDP(18); (t83 + t84) * t53 + (-MDP(13) + t80) * t52 + t77 * t50 + (-pkin(3) * MDP(15) + MDP(10) + t74) * t58 + (-MDP(11) + MDP(19) * t92 + (-t65 + t66) * MDP(20) - t75 * qJ(4)) * t57; 0; -0.2e1 * t81 + t66 * MDP(19) + MDP(12) + (t97 + t89) * pkin(3) + (0.2e1 * MDP(17) + t84 + 0.2e1 * t85 + 0.2e1 * t86) * qJ(4); t52 * MDP(18) + t75 * t58; 0; t80; MDP(18); t58 * MDP(23) + t46 * MDP(24) - t47 * MDP(25) + t79 * t57; -t77; t74; t78; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
