% Calculate joint inertia matrix for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR3_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:42:11
% EndTime: 2021-01-15 15:42:12
% DurationCPUTime: 0.18s
% Computational Cost: add. (214->63), mult. (409->101), div. (0->0), fcn. (429->6), ass. (0->33)
t60 = sin(pkin(9));
t61 = cos(pkin(9));
t63 = sin(qJ(3));
t65 = cos(qJ(3));
t53 = t60 * t65 + t61 * t63;
t71 = t53 * MDP(13);
t51 = t60 * t63 - t61 * t65;
t72 = t51 * MDP(12);
t62 = sin(qJ(5));
t64 = cos(qJ(5));
t43 = t64 * t51 + t62 * t53;
t39 = t43 * MDP(21);
t44 = -t62 * t51 + t64 * t53;
t75 = -t44 * MDP(22) - t39;
t80 = t71 + t72 - t75;
t79 = MDP(15) * pkin(3);
t59 = -t65 * pkin(3) - pkin(2);
t78 = 0.2e1 * t51 * pkin(4) + 0.2e1 * t59;
t77 = pkin(3) * t60;
t76 = -qJ(4) - pkin(6);
t58 = t61 * pkin(3) + pkin(4);
t74 = (t64 * t58 - t62 * t77) * MDP(21);
t73 = (-t62 * t58 - t64 * t77) * MDP(22);
t70 = t59 * MDP(15);
t69 = t65 * MDP(10);
t55 = t76 * t63;
t56 = t76 * t65;
t45 = t61 * t55 + t60 * t56;
t37 = -t53 * pkin(7) + t45;
t46 = t60 * t55 - t61 * t56;
t38 = -t51 * pkin(7) + t46;
t68 = t44 * MDP(18) - t43 * MDP(19) + (t64 * t37 - t62 * t38) * MDP(21) + (-t62 * t37 - t64 * t38) * MDP(22);
t1 = [MDP(1) + (t51 ^ 2 + t53 ^ 2) * MDP(15); (-t51 * t45 + t53 * t46) * MDP(15); MDP(2) + 0.2e1 * pkin(2) * t69 + 0.2e1 * (-t45 * t53 - t46 * t51) * MDP(14) + (t45 ^ 2 + t46 ^ 2) * MDP(15) + t39 * t78 + (t70 + 0.2e1 * t71 + 0.2e1 * t72) * t59 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t63 + 0.2e1 * t65 * MDP(6)) * t63 + (MDP(16) * t44 - 0.2e1 * t43 * MDP(17) + MDP(22) * t78) * t44; t69 - t63 * MDP(11) + (-t51 * t61 + t53 * t60) * t79 - t80; t45 * MDP(12) - t46 * MDP(13) + t63 * MDP(7) + t65 * MDP(8) + (-t63 * MDP(10) - t65 * MDP(11)) * pkin(6) + ((-t51 * t60 - t53 * t61) * MDP(14) + (t45 * t61 + t46 * t60) * MDP(15)) * pkin(3) + t68; MDP(9) + MDP(20) + 0.2e1 * t74 + 0.2e1 * t73 + (0.2e1 * t61 * MDP(12) - 0.2e1 * t60 * MDP(13) + (t60 ^ 2 + t61 ^ 2) * t79) * pkin(3); 0; t70 + t80; 0; MDP(15); t75; t68; MDP(20) + t73 + t74; 0; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
