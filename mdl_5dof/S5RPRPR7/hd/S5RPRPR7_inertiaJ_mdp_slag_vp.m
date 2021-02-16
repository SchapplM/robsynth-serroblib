% Calculate joint inertia matrix for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:07
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR7_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:06:27
% EndTime: 2021-01-15 12:06:28
% DurationCPUTime: 0.21s
% Computational Cost: add. (268->77), mult. (493->124), div. (0->0), fcn. (487->8), ass. (0->39)
t66 = cos(pkin(8));
t60 = -t66 * pkin(1) - pkin(2);
t87 = cos(qJ(3));
t56 = -t87 * pkin(3) + t60;
t88 = 0.2e1 * t56;
t64 = sin(pkin(8));
t80 = t64 * pkin(1) + pkin(6);
t75 = t87 * t80;
t51 = t87 * qJ(4) + t75;
t63 = sin(pkin(9));
t65 = cos(pkin(9));
t68 = sin(qJ(3));
t76 = (-qJ(4) - t80) * t68;
t47 = t63 * t51 - t65 * t76;
t55 = t63 * t87 + t65 * t68;
t86 = t47 * t55;
t67 = sin(qJ(5));
t69 = cos(qJ(5));
t85 = t67 * t69;
t84 = MDP(15) * pkin(3);
t53 = t63 * t68 - t65 * t87;
t83 = t53 * MDP(20);
t82 = t55 * MDP(13);
t81 = MDP(17) * t85;
t79 = t87 * MDP(10);
t78 = t69 * MDP(21) - t67 * MDP(22);
t77 = MDP(21) * t67 + MDP(22) * t69;
t74 = MDP(12) + t78;
t73 = (MDP(18) * t69 - MDP(19) * t67) * t55;
t72 = t67 * MDP(18) + t69 * MDP(19) - t77 * (t63 * pkin(3) + pkin(7));
t62 = t69 ^ 2;
t61 = t67 ^ 2;
t59 = -t65 * pkin(3) - pkin(4);
t52 = t55 ^ 2;
t49 = t65 * t51 + t63 * t76;
t46 = t53 * pkin(4) - t55 * pkin(7) + t56;
t45 = t67 * t46 + t69 * t49;
t44 = t69 * t46 - t67 * t49;
t1 = [MDP(1) - 0.2e1 * t60 * t79 + t82 * t88 + (t47 ^ 2 + t49 ^ 2 + t56 ^ 2) * MDP(15) + (t64 ^ 2 + t66 ^ 2) * MDP(4) * pkin(1) ^ 2 + (t62 * MDP(16) - 0.2e1 * t81) * t52 + (0.2e1 * t60 * MDP(11) + MDP(5) * t68 + 0.2e1 * t87 * MDP(6)) * t68 + (MDP(12) * t88 + 0.2e1 * t73 + t83) * t53 + 0.2e1 * (-t49 * t53 + t86) * MDP(14) + 0.2e1 * (t44 * t53 + t67 * t86) * MDP(21) + 0.2e1 * (-t45 * t53 + t69 * t86) * MDP(22); (t47 * t53 + t49 * t55) * MDP(15); MDP(4) + (t53 ^ 2 + t52) * MDP(15); t87 * MDP(8) - MDP(11) * t75 - t49 * MDP(13) + (-t80 * MDP(10) + MDP(7)) * t68 - t74 * t47 + t72 * t53 + (MDP(16) * t85 + (-t61 + t62) * MDP(17) + t77 * t59) * t55 + ((-t53 * t63 - t55 * t65) * MDP(14) + (-t47 * t65 + t49 * t63) * MDP(15)) * pkin(3); t79 - t68 * MDP(11) + (t63 * t84 - MDP(13)) * t55 + (-t65 * t84 - t74) * t53; 0.2e1 * t81 + t61 * MDP(16) + MDP(9) + (t63 ^ 2 + t65 ^ 2) * MDP(15) * pkin(3) ^ 2 - 0.2e1 * t78 * t59 + 0.2e1 * (t65 * MDP(12) - t63 * MDP(13)) * pkin(3); t56 * MDP(15) + t74 * t53 + t82; 0; 0; MDP(15); t44 * MDP(21) - t45 * MDP(22) + t73 + t83; -t77 * t55; t72; t78; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
