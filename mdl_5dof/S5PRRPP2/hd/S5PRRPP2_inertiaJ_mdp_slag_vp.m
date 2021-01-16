% Calculate joint inertia matrix for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:33
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRRPP2_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:32:34
% EndTime: 2021-01-15 15:32:35
% DurationCPUTime: 0.18s
% Computational Cost: add. (232->70), mult. (464->103), div. (0->0), fcn. (452->6), ass. (0->34)
t94 = 2 * MDP(12);
t93 = MDP(15) * pkin(3);
t92 = MDP(14) + MDP(17);
t83 = MDP(15) + MDP(19);
t65 = sin(pkin(8));
t66 = cos(pkin(8));
t67 = sin(qJ(3));
t69 = cos(qJ(3));
t55 = t65 * t67 - t66 * t69;
t56 = t65 * t69 + t66 * t67;
t84 = MDP(13) - MDP(18);
t63 = -t69 * pkin(3) - pkin(2);
t87 = t63 * MDP(15);
t45 = t55 * pkin(4) - t56 * qJ(5) + t63;
t88 = t45 * MDP(19);
t91 = (MDP(12) + MDP(16)) * t55 + t84 * t56 + t87 + t88;
t90 = 2 * MDP(16);
t89 = qJ(4) + pkin(6);
t86 = t69 * MDP(10);
t81 = t89 * t67;
t61 = t66 * pkin(3) + pkin(4);
t76 = -t61 * MDP(19) - MDP(16);
t75 = -MDP(12) + t76;
t59 = t65 * pkin(3) + qJ(5);
t74 = MDP(19) * t59 - t84;
t73 = -t67 * MDP(10) - t69 * MDP(11);
t70 = cos(qJ(2));
t68 = sin(qJ(2));
t58 = t89 * t69;
t54 = t55 * t68;
t52 = t56 * t68;
t50 = t66 * t58 - t65 * t81;
t48 = t65 * t58 + t66 * t81;
t1 = [MDP(1) + t83 * (t52 ^ 2 + t54 ^ 2 + t70 ^ 2); -t68 * MDP(4) + (-t67 * MDP(11) + MDP(3) + t86 - t91) * t70 + t83 * (t52 * t48 - t54 * t50) + t92 * (t52 * t56 + t54 * t55); MDP(2) + 0.2e1 * pkin(2) * t86 + (0.2e1 * t56 * MDP(13) + t55 * t94 + t87) * t63 + (-0.2e1 * t56 * MDP(18) + t55 * t90 + t88) * t45 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t67 + 0.2e1 * t69 * MDP(6)) * t67 + t83 * (t48 ^ 2 + t50 ^ 2) + 0.2e1 * t92 * (t48 * t56 - t50 * t55); t73 * t68 - t74 * t54 + t75 * t52 + (-t52 * t66 - t54 * t65) * t93; t67 * MDP(7) + t69 * MDP(8) + (-t59 * t55 - t61 * t56) * MDP(17) + t74 * t50 + t75 * t48 + t73 * pkin(6) + ((-t55 * t65 - t56 * t66) * MDP(14) + (-t48 * t66 + t50 * t65) * MDP(15)) * pkin(3); MDP(9) + (t59 ^ 2 + t61 ^ 2) * MDP(19) + t61 * t90 + 0.2e1 * t59 * MDP(18) + (t66 * t94 - 0.2e1 * t65 * MDP(13) + (t65 ^ 2 + t66 ^ 2) * t93) * pkin(3); -t83 * t70; t91; 0; t83; t52 * MDP(19); t56 * MDP(17) + t48 * MDP(19); t76; 0; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
