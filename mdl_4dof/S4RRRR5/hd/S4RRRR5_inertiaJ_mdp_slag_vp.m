% Calculate joint inertia matrix for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S4RRRR5_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:14
% EndTime: 2019-12-31 17:28:15
% DurationCPUTime: 0.29s
% Computational Cost: add. (215->84), mult. (453->125), div. (0->0), fcn. (417->6), ass. (0->45)
t79 = sin(qJ(3));
t82 = cos(qJ(3));
t78 = sin(qJ(4));
t81 = cos(qJ(4));
t66 = t78 * t79 - t81 * t82;
t67 = t78 * t82 + t79 * t81;
t102 = pkin(6) + pkin(7);
t70 = t102 * t79;
t71 = t102 * t82;
t89 = t67 * MDP(20) - t66 * MDP(21) + (-t70 * t81 - t71 * t78) * MDP(23) - (-t70 * t78 + t71 * t81) * MDP(24);
t108 = -(MDP(16) * t79 + MDP(17) * t82) * pkin(6) + t79 * MDP(13) + t82 * MDP(14) + t89;
t84 = (MDP(23) * t81 - MDP(24) * t78) * pkin(3);
t80 = sin(qJ(2));
t61 = t67 * t80;
t62 = t66 * t80;
t96 = -t62 * MDP(20) - t61 * MDP(21);
t101 = pkin(5) * t79;
t83 = cos(qJ(2));
t69 = -pkin(2) * t83 - t80 * pkin(6) - pkin(1);
t65 = t82 * t69;
t98 = pkin(7) * t80;
t51 = -t82 * t98 + t65 + (-pkin(3) - t101) * t83;
t99 = pkin(5) * t83;
t91 = t82 * t99;
t52 = t91 + (t69 - t98) * t79;
t48 = t81 * t51 - t78 * t52;
t49 = t78 * t51 + t81 * t52;
t107 = t48 * MDP(23) - t49 * MDP(24) + t96;
t105 = -2 * MDP(19);
t104 = 0.2e1 * MDP(23);
t103 = 0.2e1 * MDP(24);
t100 = pkin(5) * t82;
t97 = t79 * t82;
t95 = MDP(18) * t67;
t92 = MDP(15) + MDP(22);
t90 = MDP(12) * t97;
t88 = MDP(13) * t82 - MDP(14) * t79;
t76 = t82 ^ 2;
t75 = t80 ^ 2;
t74 = t79 ^ 2;
t73 = -pkin(3) * t82 - pkin(2);
t68 = (pkin(3) * t79 + pkin(5)) * t80;
t58 = t79 * t69 + t91;
t57 = -t79 * t99 + t65;
t1 = [-0.2e1 * pkin(1) * t80 * MDP(10) + MDP(1) + t92 * t83 ^ 2 - (-MDP(18) * t62 + t61 * t105) * t62 + (MDP(11) * t76 + MDP(4) - 0.2e1 * t90) * t75 + 0.2e1 * (pkin(1) * MDP(9) + (MDP(5) - t88) * t80 - t96) * t83 + 0.2e1 * (t75 * t101 - t57 * t83) * MDP(16) + 0.2e1 * (t75 * t100 + t58 * t83) * MDP(17) + (-t48 * t83 + t68 * t61) * t104 + (t49 * t83 - t68 * t62) * t103; -t62 * t95 + (-t61 * t67 + t62 * t66) * MDP(19) + (t73 * t61 + t68 * t66) * MDP(23) + (-t73 * t62 + t68 * t67) * MDP(24) + (-pkin(5) * MDP(10) + MDP(7) - t108) * t83 + (MDP(6) - pkin(5) * MDP(9) + MDP(11) * t97 + (-t74 + t76) * MDP(12) + (-pkin(2) * t79 - t100) * MDP(16) + (-pkin(2) * t82 + t101) * MDP(17)) * t80; 0.2e1 * t90 + t66 * t73 * t104 + MDP(11) * t74 + MDP(8) + 0.2e1 * (MDP(16) * t82 - MDP(17) * t79) * pkin(2) + (t73 * t103 + t66 * t105 + t95) * t67; t57 * MDP(16) - t58 * MDP(17) + (-t92 - t84) * t83 + t88 * t80 + t107; t108; 0.2e1 * t84 + t92; -t83 * MDP(22) + t107; t89; MDP(22) + t84; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
