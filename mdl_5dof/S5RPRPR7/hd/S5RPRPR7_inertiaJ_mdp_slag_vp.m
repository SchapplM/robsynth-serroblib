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
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRPR7_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:51
% EndTime: 2019-12-31 18:19:52
% DurationCPUTime: 0.18s
% Computational Cost: add. (238->68), mult. (443->114), div. (0->0), fcn. (441->8), ass. (0->36)
t61 = sin(pkin(9));
t63 = cos(pkin(9));
t66 = sin(qJ(3));
t80 = cos(qJ(3));
t50 = t61 * t66 - t63 * t80;
t65 = sin(qJ(5));
t67 = cos(qJ(5));
t73 = MDP(19) * t67 - t65 * MDP(20);
t82 = t73 * t50;
t62 = sin(pkin(8));
t56 = pkin(1) * t62 + pkin(6);
t81 = qJ(4) + t56;
t48 = t81 * t80;
t75 = t81 * t66;
t44 = t48 * t61 + t63 * t75;
t52 = t61 * t80 + t63 * t66;
t79 = t44 * t52;
t78 = t65 * t67;
t77 = MDP(18) * t50;
t76 = MDP(15) * t78;
t64 = cos(pkin(8));
t58 = -t64 * pkin(1) - pkin(2);
t74 = t80 * MDP(10);
t72 = MDP(19) * t65 + MDP(20) * t67;
t53 = -t80 * pkin(3) + t58;
t71 = (MDP(16) * t67 - MDP(17) * t65) * t52;
t70 = MDP(16) * t65 + MDP(17) * t67 - t72 * (pkin(3) * t61 + pkin(7));
t60 = t67 ^ 2;
t59 = t65 ^ 2;
t57 = -pkin(3) * t63 - pkin(4);
t49 = t52 ^ 2;
t46 = t63 * t48 - t61 * t75;
t43 = t50 * pkin(4) - t52 * pkin(7) + t53;
t42 = t43 * t65 + t46 * t67;
t41 = t67 * t43 - t65 * t46;
t1 = [MDP(1) - 0.2e1 * t58 * t74 + (t44 ^ 2 + t46 ^ 2 + t53 ^ 2) * MDP(13) + (t62 ^ 2 + t64 ^ 2) * MDP(4) * pkin(1) ^ 2 + (t60 * MDP(14) - 0.2e1 * t76) * t49 + (0.2e1 * t58 * MDP(11) + MDP(5) * t66 + 0.2e1 * t80 * MDP(6)) * t66 + (0.2e1 * t71 + t77) * t50 + 0.2e1 * (-t46 * t50 + t79) * MDP(12) + 0.2e1 * (t41 * t50 + t65 * t79) * MDP(19) + 0.2e1 * (-t42 * t50 + t67 * t79) * MDP(20); (t44 * t50 + t46 * t52) * MDP(13); MDP(4) + (t50 ^ 2 + t49) * MDP(13); t66 * MDP(7) + t80 * MDP(8) + (-t66 * MDP(10) - t80 * MDP(11)) * t56 - t73 * t44 + t70 * t50 + (MDP(14) * t78 + (-t59 + t60) * MDP(15) + t72 * t57) * t52 + ((-t50 * t61 - t52 * t63) * MDP(12) + (-t44 * t63 + t46 * t61) * MDP(13)) * pkin(3); t74 - t66 * MDP(11) - t82 + (-t50 * t63 + t52 * t61) * MDP(13) * pkin(3); 0.2e1 * t76 + t59 * MDP(14) + MDP(9) + (t61 ^ 2 + t63 ^ 2) * MDP(13) * pkin(3) ^ 2 - 0.2e1 * t73 * t57; t53 * MDP(13) + t82; 0; 0; MDP(13); MDP(19) * t41 - MDP(20) * t42 + t71 + t77; -t72 * t52; t70; t73; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
