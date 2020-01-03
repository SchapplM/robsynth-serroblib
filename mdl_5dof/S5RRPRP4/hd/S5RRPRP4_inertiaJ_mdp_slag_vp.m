% Calculate joint inertia matrix for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RRPRP4_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:56
% EndTime: 2019-12-31 19:52:57
% DurationCPUTime: 0.20s
% Computational Cost: add. (193->76), mult. (270->91), div. (0->0), fcn. (144->4), ass. (0->39)
t57 = sin(qJ(4));
t59 = cos(qJ(4));
t80 = 2 * MDP(8);
t84 = 0.2e1 * t57 * MDP(15) + 0.2e1 * t59 * MDP(16) + t80;
t60 = cos(qJ(2));
t49 = -t60 * pkin(1) - pkin(2);
t83 = t49 * MDP(9);
t55 = t59 ^ 2;
t44 = t57 ^ 2 + t55;
t68 = t44 * MDP(18);
t67 = t44 * MDP(20);
t82 = -0.2e1 * t59;
t81 = -2 * MDP(7);
t79 = 2 * MDP(17);
t78 = -0.2e1 * MDP(18);
t58 = sin(qJ(2));
t77 = t58 * pkin(1);
t76 = pkin(2) * MDP(9);
t40 = t57 * pkin(4) - t59 * qJ(5) + qJ(3);
t37 = t40 + t77;
t75 = t37 + t40;
t74 = MDP(7) - t68;
t47 = qJ(3) + t77;
t73 = qJ(3) + t47;
t72 = t37 * MDP(20);
t71 = t59 * MDP(20);
t70 = -MDP(16) + MDP(19);
t69 = t57 * MDP(11) * t82 + t55 * MDP(10) + MDP(4);
t66 = -pkin(4) * MDP(20) - MDP(17);
t41 = t59 * pkin(4) + t57 * qJ(5);
t65 = t59 * MDP(12) - t57 * MDP(13) - t41 * MDP(18);
t63 = MDP(19) * t82 + t57 * t79;
t62 = (MDP(15) - t66) * t59 + (MDP(20) * qJ(5) + t70) * t57;
t61 = -pkin(2) - pkin(7);
t51 = t59 * MDP(18);
t46 = -pkin(7) + t49;
t38 = t44 * t61;
t36 = t44 * t46;
t1 = [MDP(1) + t46 ^ 2 * t67 + (t63 + t72) * t37 + 0.2e1 * (t60 * MDP(5) - t58 * MDP(6)) * pkin(1) + t36 * t78 + t69 + ((2 * MDP(7)) + t83) * t49 + (t47 * MDP(9) + t84) * t47; pkin(2) * t81 + qJ(3) * t80 + (-t49 * pkin(2) + t47 * qJ(3)) * MDP(9) + t40 * t72 - t61 * t68 + (t61 * t67 - t68) * t46 + ((MDP(5) - MDP(7)) * t60 + (-MDP(6) + MDP(8)) * t58) * pkin(1) + (t73 * MDP(16) - t75 * MDP(19)) * t59 + (t73 * MDP(15) + t75 * MDP(17)) * t57 + t69; t38 * t78 + t61 ^ 2 * t67 + (MDP(20) * t40 + t63) * t40 + (t81 + t76) * pkin(2) + (MDP(9) * qJ(3) + t84) * qJ(3) + t69; t36 * MDP(20) + t74 + t83; t38 * MDP(20) + t74 - t76; MDP(9) + t67; t46 * t62 + t65; t61 * t62 + t65; t41 * MDP(20) + (MDP(15) + MDP(17)) * t59 + t70 * t57; MDP(14) + pkin(4) * t79 + 0.2e1 * qJ(5) * MDP(19) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(20); -t46 * t71 + t51; -t61 * t71 + t51; -t71; t66; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
