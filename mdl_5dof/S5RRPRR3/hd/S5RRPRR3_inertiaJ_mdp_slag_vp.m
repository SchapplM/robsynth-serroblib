% Calculate joint inertia matrix for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRPRR3_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:35
% EndTime: 2020-01-03 12:00:36
% DurationCPUTime: 0.14s
% Computational Cost: add. (202->61), mult. (341->74), div. (0->0), fcn. (270->8), ass. (0->40)
t63 = sin(qJ(5));
t66 = cos(qJ(5));
t78 = t66 * MDP(16);
t73 = -t63 * MDP(17) + t78;
t70 = 0.2e1 * t73;
t65 = sin(qJ(2));
t87 = pkin(1) * t65;
t61 = sin(pkin(9));
t86 = pkin(2) * t61;
t85 = pkin(4) * t63;
t62 = cos(pkin(9));
t55 = t62 * pkin(2) + pkin(3);
t64 = sin(qJ(4));
t84 = t64 * t55;
t83 = t63 * MDP(13) + t66 * MDP(14);
t68 = cos(qJ(2));
t56 = t68 * pkin(1) + pkin(2);
t48 = t62 * t56 - t61 * t87;
t45 = pkin(3) + t48;
t50 = t61 * t56 + t62 * t87;
t67 = cos(qJ(4));
t76 = -t67 * t45 + t64 * t50;
t82 = t76 * MDP(9);
t53 = t67 * t55;
t49 = -t64 * t86 + t53;
t81 = t49 * MDP(9);
t74 = t64 * t45 + t67 * t50;
t80 = t74 * MDP(10);
t51 = -t67 * t86 - t84;
t79 = t51 * MDP(10);
t77 = MDP(8) + (MDP(11) * t63 + 0.2e1 * MDP(12) * t66) * t63;
t75 = MDP(4) + t77;
t72 = -MDP(16) * t63 - MDP(17) * t66;
t71 = (t68 * MDP(5) - t65 * MDP(6)) * pkin(1);
t60 = pkin(4) * t66;
t46 = -pkin(4) - t49;
t44 = t46 * t63;
t39 = -pkin(4) + t76;
t38 = t39 * t63;
t1 = [MDP(1) + (t48 ^ 2 + t50 ^ 2) * MDP(7) - t39 * t70 + 0.2e1 * t71 - 0.2e1 * t80 - 0.2e1 * t82 + t75; (t53 - t76) * MDP(9) + (-t74 - t84) * MDP(10) + (t44 + t38) * MDP(17) + (-t39 - t46) * t78 + (t48 * t62 * MDP(7) + (-t67 * MDP(10) + t50 * MDP(7) - t64 * MDP(9)) * t61) * pkin(2) + t71 + t75; (t61 ^ 2 + t62 ^ 2) * MDP(7) * pkin(2) ^ 2 - t46 * t70 + 0.2e1 * t79 + 0.2e1 * t81 + t75; 0; 0; MDP(7); -t82 - t80 + (-t39 * t66 + t60) * MDP(16) + (t38 - t85) * MDP(17) + t77; t81 + t79 + (-t46 * t66 + t60) * MDP(16) + (t44 - t85) * MDP(17) + t77; 0; pkin(4) * t70 + t77; t72 * (pkin(8) + t74) + t83; t72 * (pkin(8) - t51) + t83; t73; pkin(8) * t72 + t83; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
