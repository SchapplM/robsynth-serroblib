% Calculate joint inertia matrix for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:14
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPRPP1_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:14:04
% EndTime: 2021-01-15 11:14:05
% DurationCPUTime: 0.20s
% Computational Cost: add. (265->67), mult. (445->100), div. (0->0), fcn. (422->6), ass. (0->31)
t87 = 2 * MDP(12);
t86 = MDP(15) * pkin(3);
t85 = MDP(15) + MDP(19);
t84 = 2 * MDP(16);
t83 = cos(qJ(3));
t61 = sin(pkin(8));
t63 = cos(pkin(8));
t65 = sin(qJ(3));
t52 = t61 * t65 - t63 * t83;
t54 = t61 * t83 + t63 * t65;
t64 = cos(pkin(7));
t60 = -t64 * pkin(1) - pkin(2);
t55 = -pkin(3) * t83 + t60;
t43 = t52 * pkin(4) - t54 * qJ(5) + t55;
t82 = t43 * MDP(19);
t81 = t55 * MDP(15);
t80 = MDP(13) - MDP(18);
t62 = sin(pkin(7));
t77 = t62 * pkin(1) + pkin(6);
t73 = t83 * MDP(10);
t58 = t63 * pkin(3) + pkin(4);
t72 = -t58 * MDP(19) - MDP(16);
t71 = -MDP(12) + t72;
t56 = t61 * pkin(3) + qJ(5);
t70 = MDP(19) * t56 - t80;
t69 = (-qJ(4) - t77) * t65;
t68 = t83 * t77;
t50 = qJ(4) * t83 + t68;
t47 = t63 * t50 + t61 * t69;
t45 = t61 * t50 - t63 * t69;
t1 = [MDP(1) - 0.2e1 * t60 * t73 + (t62 ^ 2 + t64 ^ 2) * MDP(4) * pkin(1) ^ 2 + (0.2e1 * t54 * MDP(13) + t52 * t87 + t81) * t55 + (-0.2e1 * t54 * MDP(18) + t52 * t84 + t82) * t43 + (0.2e1 * t60 * MDP(11) + MDP(5) * t65 + 0.2e1 * MDP(6) * t83) * t65 + t85 * (t45 ^ 2 + t47 ^ 2) + 0.2e1 * (MDP(14) + MDP(17)) * (t45 * t54 - t47 * t52); t85 * (t45 * t52 + t47 * t54); MDP(4) + t85 * (t52 ^ 2 + t54 ^ 2); t83 * MDP(8) - MDP(11) * t68 + (-t56 * t52 - t58 * t54) * MDP(17) + (-MDP(10) * t77 + MDP(7)) * t65 + t70 * t47 + t71 * t45 + ((-t52 * t61 - t54 * t63) * MDP(14) + (-t45 * t63 + t47 * t61) * MDP(15)) * pkin(3); t73 - t65 * MDP(11) + t70 * t54 + t71 * t52 + (-t52 * t63 + t54 * t61) * t86; MDP(9) + (t56 ^ 2 + t58 ^ 2) * MDP(19) + t58 * t84 + 0.2e1 * t56 * MDP(18) + (t63 * t87 - 0.2e1 * t61 * MDP(13) + (t61 ^ 2 + t63 ^ 2) * t86) * pkin(3); t81 + t82 + t80 * t54 + (MDP(12) + MDP(16)) * t52; 0; 0; t85; t54 * MDP(17) + t45 * MDP(19); t52 * MDP(19); t72; 0; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
