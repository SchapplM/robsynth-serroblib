% Calculate joint inertia matrix for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RPPRRP5_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:08:48
% EndTime: 2019-03-09 02:08:48
% DurationCPUTime: 0.27s
% Computational Cost: add. (248->109), mult. (395->146), div. (0->0), fcn. (298->4), ass. (0->48)
t70 = cos(qJ(4));
t97 = -t70 * MDP(16) - MDP(8);
t95 = 2 * MDP(24);
t64 = -pkin(7) + qJ(2);
t67 = sin(qJ(5));
t94 = t64 * t67;
t69 = cos(qJ(5));
t93 = t64 * t69;
t92 = t67 * t69;
t65 = pkin(1) + qJ(3);
t91 = -qJ(6) - pkin(8);
t60 = t67 ^ 2;
t62 = t69 ^ 2;
t90 = t60 + t62;
t68 = sin(qJ(4));
t61 = t68 ^ 2;
t63 = t70 ^ 2;
t89 = -t61 - t63;
t88 = MDP(25) * pkin(5);
t87 = qJ(6) * t70;
t86 = MDP(23) * t69;
t55 = pkin(4) * t68 - pkin(8) * t70 + t65;
t53 = t69 * t55;
t49 = -t69 * t87 + t53 + (pkin(5) - t94) * t68;
t85 = t49 * MDP(25);
t59 = -pkin(5) * t69 - pkin(4);
t84 = t59 * MDP(25);
t83 = t67 * MDP(20);
t82 = t67 * MDP(23);
t81 = t70 * MDP(25);
t80 = t68 * t93;
t79 = MDP(18) * t92;
t78 = -MDP(24) * pkin(5) + MDP(19);
t77 = -MDP(22) - t88;
t50 = t80 + (t55 - t87) * t67;
t76 = -t49 * t69 - t50 * t67;
t56 = t91 * t67;
t57 = t91 * t69;
t75 = -t56 * t69 + t57 * t67;
t74 = -t56 * t67 - t57 * t69;
t73 = t69 * MDP(22) - t82;
t72 = -MDP(15) - t73;
t71 = (qJ(2) ^ 2);
t58 = t62 * t70;
t54 = (pkin(5) * t67 - t64) * t70;
t52 = t55 * t67 + t80;
t51 = -t68 * t94 + t53;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2 + t71) * MDP(6)) + (t65 ^ 2 + t71) * MDP(9) + t61 * MDP(21) + (t49 ^ 2 + t50 ^ 2 + t54 ^ 2) * MDP(25) + (t62 * MDP(17) + MDP(10) - 0.2e1 * t79) * t63 + 0.2e1 * (t51 * t68 - t63 * t94) * MDP(22) + 0.2e1 * (-t52 * t68 - t63 * t93) * MDP(23) + (2 * (MDP(5) + MDP(7)) * qJ(2)) + 0.2e1 * (t68 * MDP(15) - t97) * t65 + (t76 * t95 + 0.2e1 * (t69 * MDP(19) - MDP(11) - t83) * t68) * t70; MDP(4) - (pkin(1) * MDP(6)) - t65 * MDP(9) + (t60 * t70 + t58) * MDP(24) + t76 * MDP(25) + t72 * t68 + t97; t90 * MDP(25) + MDP(6) + MDP(9); -t54 * t81 + qJ(2) * MDP(9) + MDP(7) + (t50 * t68 * MDP(25) + t89 * MDP(23)) * t69 + (t89 * MDP(22) - t68 * t85) * t67; 0; MDP(9) + (t90 * t61 + t63) * MDP(25); t58 * MDP(18) + (-t49 * t67 + t50 * t69) * MDP(24) + (t49 * t56 - t50 * t57 + t54 * t59) * MDP(25) + (-t64 * MDP(16) + t67 * MDP(19) + t69 * MDP(20) - MDP(13) + (-MDP(22) * t67 - t86) * pkin(8)) * t68 + (MDP(12) + t64 * MDP(15) + MDP(17) * t92 - t60 * MDP(18) + (-pkin(4) * t67 + t93) * MDP(22) + (-pkin(4) * t69 - t94) * MDP(23) + t75 * MDP(24)) * t70; t75 * MDP(25); (-t72 - t84) * t70 + (t90 * MDP(24) + t74 * MDP(25) - MDP(16)) * t68; MDP(14) + t60 * MDP(17) + 0.2e1 * t79 + t74 * t95 + (t56 ^ 2 + t57 ^ 2 + t59 ^ 2) * MDP(25) + 0.2e1 * t73 * pkin(4); pkin(5) * t85 + MDP(21) * t68 + t51 * MDP(22) - t52 * MDP(23) + (t78 * t69 - t83) * t70; t77 * t69 + t82; (t77 * t67 - t86) * t68; t56 * t88 + (-MDP(23) * pkin(8) + MDP(20)) * t69 + (-MDP(22) * pkin(8) + t78) * t67; MDP(25) * pkin(5) ^ 2 + MDP(21); t54 * MDP(25); 0; -t81; t84; 0; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
