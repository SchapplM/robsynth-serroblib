% Calculate joint inertia matrix for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP5_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:58:29
% EndTime: 2019-12-31 20:58:30
% DurationCPUTime: 0.25s
% Computational Cost: add. (326->105), mult. (534->136), div. (0->0), fcn. (468->4), ass. (0->37)
t93 = 2 * qJ(4);
t77 = cos(qJ(2));
t92 = 0.2e1 * t77;
t91 = 2 * MDP(18);
t90 = 2 * MDP(22);
t78 = pkin(3) + pkin(4);
t89 = -pkin(7) - pkin(6);
t74 = sin(qJ(3));
t75 = sin(qJ(2));
t76 = cos(qJ(3));
t58 = t74 * t75 - t76 * t77;
t73 = t74 * pkin(2);
t67 = t73 + qJ(4);
t88 = t67 * t58;
t87 = qJ(4) * t58;
t86 = t74 * MDP(17);
t85 = MDP(18) + MDP(22);
t84 = MDP(20) + MDP(23);
t72 = -t77 * pkin(2) - pkin(1);
t70 = t76 * pkin(2) + pkin(3);
t61 = t89 * t75;
t62 = t89 * t77;
t54 = -t76 * t61 - t74 * t62;
t55 = t74 * t61 - t76 * t62;
t59 = t74 * t77 + t76 * t75;
t83 = t59 * qJ(4) - t72;
t47 = -t59 * qJ(5) + t54;
t48 = t58 * qJ(5) + t55;
t82 = t59 * MDP(13) - t58 * MDP(14) - t47 * MDP(22) + t48 * MDP(23) + (-MDP(17) + MDP(20)) * t55 + (-MDP(16) - MDP(18)) * t54;
t81 = 0.2e1 * pkin(3);
t80 = qJ(4) ^ 2;
t65 = pkin(4) + t70;
t64 = t67 ^ 2;
t63 = t67 * qJ(4);
t49 = t58 * pkin(3) - t83;
t44 = -t78 * t58 + t83;
t1 = [MDP(1) + pkin(1) * MDP(9) * t92 + (t49 ^ 2 + t54 ^ 2 + t55 ^ 2) * MDP(21) + (t44 ^ 2 + t47 ^ 2 + t48 ^ 2) * MDP(25) + 0.2e1 * (t72 * MDP(16) + t49 * MDP(18) - t44 * MDP(22)) * t58 + (MDP(11) * t59 - 0.2e1 * t58 * MDP(12) + 0.2e1 * t72 * MDP(17) - 0.2e1 * t49 * MDP(20) + 0.2e1 * t44 * MDP(23)) * t59 + 0.2e1 * (t54 * t59 - t55 * t58) * MDP(19) + 0.2e1 * (-t47 * t59 + t48 * t58) * MDP(24) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t75 + MDP(5) * t92) * t75; t75 * MDP(6) + t77 * MDP(7) + (-t70 * t59 - t88) * MDP(19) + (-t54 * t70 + t55 * t67) * MDP(21) + (t65 * t59 + t88) * MDP(24) + (-t47 * t65 + t48 * t67) * MDP(25) + (-t77 * MDP(10) - t75 * MDP(9)) * pkin(6) + t82; MDP(8) + MDP(15) + (t70 ^ 2 + t64) * MDP(21) + (t65 ^ 2 + t64) * MDP(25) + 0.2e1 * (t76 * MDP(16) - t86) * pkin(2) + t70 * t91 + t65 * t90 + 0.2e1 * t84 * t67; (-pkin(3) * t59 - t87) * MDP(19) + (-t54 * pkin(3) + t55 * qJ(4)) * MDP(21) + (t78 * t59 + t87) * MDP(24) + (t48 * qJ(4) - t47 * t78) * MDP(25) + t82; MDP(15) + t81 * MDP(18) + (t70 * pkin(3) + t63) * MDP(21) + (0.2e1 * pkin(4) + t81) * MDP(22) + (t65 * t78 + t63) * MDP(25) + t84 * (t93 + t73) + (-t86 + (MDP(16) + t85) * t76) * pkin(2); MDP(15) + pkin(3) * t91 + (pkin(3) ^ 2 + t80) * MDP(21) + t78 * t90 + (t78 ^ 2 + t80) * MDP(25) + t84 * t93; t54 * MDP(21) + t47 * MDP(25) + (MDP(19) - MDP(24)) * t59; -t70 * MDP(21) - t65 * MDP(25) - t85; -pkin(3) * MDP(21) - t78 * MDP(25) - t85; MDP(21) + MDP(25); -t58 * MDP(22) + t59 * MDP(23) + t44 * MDP(25); 0; 0; 0; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
