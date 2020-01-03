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
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRP5_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:57
% EndTime: 2019-12-31 19:54:58
% DurationCPUTime: 0.26s
% Computational Cost: add. (448->89), mult. (802->129), div. (0->0), fcn. (841->6), ass. (0->37)
t94 = MDP(18) + MDP(20);
t93 = MDP(19) - MDP(22);
t92 = 2 * pkin(4);
t80 = cos(qJ(2));
t91 = 0.2e1 * t80;
t90 = 2 * MDP(22);
t76 = sin(pkin(8));
t89 = pkin(2) * t76;
t88 = cos(qJ(4));
t87 = -qJ(3) - pkin(6);
t79 = sin(qJ(2));
t69 = t87 * t79;
t70 = t87 * t80;
t77 = cos(pkin(8));
t57 = t76 * t69 - t77 * t70;
t74 = t77 * pkin(2) + pkin(3);
t78 = sin(qJ(4));
t86 = t78 * t74 + t88 * t89;
t61 = t88 * t74 - t78 * t89;
t85 = t61 * MDP(18);
t84 = t86 * MDP(19);
t75 = -t80 * pkin(2) - pkin(1);
t56 = t77 * t69 + t76 * t70;
t64 = -t76 * t79 + t77 * t80;
t58 = -t64 * pkin(3) + t75;
t51 = t64 * pkin(7) + t57;
t65 = t76 * t80 + t77 * t79;
t82 = -t65 * pkin(7) + t56;
t48 = t78 * t51 - t88 * t82;
t49 = t88 * t51 + t78 * t82;
t54 = -t88 * t64 + t78 * t65;
t55 = t78 * t64 + t88 * t65;
t83 = t55 * MDP(15) - t54 * MDP(16) - t48 * t94 - t93 * t49;
t60 = -pkin(4) - t61;
t59 = qJ(5) + t86;
t47 = t54 * pkin(4) - t55 * qJ(5) + t58;
t1 = [MDP(1) + pkin(1) * MDP(9) * t91 + 0.2e1 * (-t56 * t65 + t57 * t64) * MDP(11) + (t56 ^ 2 + t57 ^ 2 + t75 ^ 2) * MDP(12) + (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) * MDP(23) + 0.2e1 * (t58 * MDP(18) + t47 * MDP(20) - MDP(21) * t49) * t54 + (MDP(13) * t55 - 0.2e1 * t54 * MDP(14) + 0.2e1 * t58 * MDP(19) + 0.2e1 * MDP(21) * t48 - 0.2e1 * t47 * MDP(22)) * t55 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t79 + MDP(5) * t91) * t79; t79 * MDP(6) + t80 * MDP(7) + (-t59 * t54 + t60 * t55) * MDP(21) + (t48 * t60 + t49 * t59) * MDP(23) + (-t80 * MDP(10) - t79 * MDP(9)) * pkin(6) + ((t64 * t76 - t65 * t77) * MDP(11) + (t56 * t77 + t57 * t76) * MDP(12)) * pkin(2) + t83; MDP(8) + MDP(17) + (t59 ^ 2 + t60 ^ 2) * MDP(23) + (t76 ^ 2 + t77 ^ 2) * MDP(12) * pkin(2) ^ 2 + 0.2e1 * t85 - 0.2e1 * t84 - 0.2e1 * t60 * MDP(20) + t59 * t90; t75 * MDP(12) + t47 * MDP(23) + t54 * t94 + t93 * t55; 0; MDP(12) + MDP(23); (-pkin(4) * t55 - t54 * qJ(5)) * MDP(21) + (-t48 * pkin(4) + t49 * qJ(5)) * MDP(23) + t83; MDP(17) + t85 - t84 + (t92 + t61) * MDP(20) + (0.2e1 * qJ(5) + t86) * MDP(22) + (-t60 * pkin(4) + t59 * qJ(5)) * MDP(23); 0; MDP(17) + MDP(20) * t92 + qJ(5) * t90 + ((pkin(4) ^ 2) + qJ(5) ^ 2) * MDP(23); t55 * MDP(21) + t48 * MDP(23); t60 * MDP(23) - MDP(20); 0; -MDP(23) * pkin(4) - MDP(20); MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
