% Calculate joint inertia matrix for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:19
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP5_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:18:21
% EndTime: 2021-01-15 20:18:22
% DurationCPUTime: 0.28s
% Computational Cost: add. (468->97), mult. (840->139), div. (0->0), fcn. (873->6), ass. (0->40)
t96 = MDP(20) + MDP(22);
t95 = MDP(21) - MDP(24);
t94 = 2 * pkin(4);
t79 = cos(qJ(2));
t93 = 0.2e1 * t79;
t92 = 2 * MDP(24);
t75 = sin(pkin(8));
t91 = pkin(2) * t75;
t90 = cos(qJ(4));
t89 = -qJ(3) - pkin(6);
t76 = cos(pkin(8));
t73 = t76 * pkin(2) + pkin(3);
t77 = sin(qJ(4));
t88 = t77 * t73 + t90 * t91;
t61 = t90 * t73 - t77 * t91;
t87 = t61 * MDP(20);
t86 = t88 * MDP(21);
t78 = sin(qJ(2));
t64 = t75 * t78 - t76 * t79;
t85 = t64 * MDP(11);
t65 = t75 * t79 + t76 * t78;
t84 = t65 * MDP(12);
t74 = -t79 * pkin(2) - pkin(1);
t83 = t74 * MDP(14);
t67 = t89 * t78;
t68 = t89 * t79;
t56 = t76 * t67 + t75 * t68;
t57 = t75 * t67 - t76 * t68;
t58 = t64 * pkin(3) + t74;
t51 = -t64 * pkin(7) + t57;
t81 = -t65 * pkin(7) + t56;
t48 = t77 * t51 - t90 * t81;
t49 = t90 * t51 + t77 * t81;
t54 = t90 * t64 + t77 * t65;
t55 = -t77 * t64 + t90 * t65;
t82 = t55 * MDP(17) - t54 * MDP(18) - t96 * t48 - t95 * t49;
t60 = -pkin(4) - t61;
t59 = qJ(5) + t88;
t47 = t54 * pkin(4) - t55 * qJ(5) + t58;
t1 = [MDP(1) + pkin(1) * MDP(9) * t93 + 0.2e1 * (-t56 * t65 - t57 * t64) * MDP(13) + (t56 ^ 2 + t57 ^ 2) * MDP(14) + (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) * MDP(25) + (t83 + 0.2e1 * t84 + 0.2e1 * t85) * t74 + 0.2e1 * (t58 * MDP(20) + t47 * MDP(22) - MDP(23) * t49) * t54 + (MDP(15) * t55 - 0.2e1 * t54 * MDP(16) + 0.2e1 * t58 * MDP(21) + 0.2e1 * MDP(23) * t48 - 0.2e1 * t47 * MDP(24)) * t55 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t78 + MDP(5) * t93) * t78; t78 * MDP(6) + t79 * MDP(7) + t56 * MDP(11) - t57 * MDP(12) + (-t59 * t54 + t60 * t55) * MDP(23) + (t48 * t60 + t49 * t59) * MDP(25) + (-t79 * MDP(10) - t78 * MDP(9)) * pkin(6) + ((-t64 * t75 - t65 * t76) * MDP(13) + (t56 * t76 + t57 * t75) * MDP(14)) * pkin(2) + t82; MDP(8) + MDP(19) + (t59 ^ 2 + t60 ^ 2) * MDP(25) + 0.2e1 * t87 - 0.2e1 * t86 - 0.2e1 * t60 * MDP(22) + t59 * t92 + (0.2e1 * t76 * MDP(11) - 0.2e1 * t75 * MDP(12) + (t75 ^ 2 + t76 ^ 2) * MDP(14) * pkin(2)) * pkin(2); t47 * MDP(25) + t96 * t54 + t95 * t55 + t83 + t84 + t85; 0; MDP(14) + MDP(25); (-pkin(4) * t55 - t54 * qJ(5)) * MDP(23) + (-t48 * pkin(4) + t49 * qJ(5)) * MDP(25) + t82; MDP(19) + t87 - t86 + (t94 + t61) * MDP(22) + (0.2e1 * qJ(5) + t88) * MDP(24) + (-t60 * pkin(4) + t59 * qJ(5)) * MDP(25); 0; MDP(19) + MDP(22) * t94 + qJ(5) * t92 + ((pkin(4) ^ 2) + qJ(5) ^ 2) * MDP(25); t55 * MDP(23) + t48 * MDP(25); t60 * MDP(25) - MDP(22); 0; -MDP(25) * pkin(4) - MDP(22); MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
