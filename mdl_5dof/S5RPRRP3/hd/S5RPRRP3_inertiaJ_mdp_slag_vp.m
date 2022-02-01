% Calculate joint inertia matrix for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP3_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:30:20
% EndTime: 2022-01-23 09:30:20
% DurationCPUTime: 0.18s
% Computational Cost: add. (248->70), mult. (417->101), div. (0->0), fcn. (398->6), ass. (0->36)
t60 = sin(pkin(8));
t55 = t60 * pkin(1) + pkin(6);
t81 = pkin(7) + t55;
t62 = sin(qJ(4));
t63 = sin(qJ(3));
t75 = cos(qJ(4));
t76 = cos(qJ(3));
t51 = t62 * t63 - t75 * t76;
t42 = t51 * MDP(19);
t53 = t62 * t76 + t63 * t75;
t46 = t53 * MDP(20);
t80 = t42 + t46;
t79 = t51 * MDP(17) + t53 * MDP(18);
t78 = 0.2e1 * MDP(19);
t77 = pkin(3) * t62;
t74 = MDP(22) * pkin(4);
t61 = cos(pkin(8));
t56 = -t61 * pkin(1) - pkin(2);
t65 = -pkin(3) * t76 + t56;
t40 = t51 * pkin(4) + t65;
t73 = t40 * MDP(22);
t59 = t75 * pkin(3);
t58 = t59 + pkin(4);
t72 = t58 * MDP(22);
t48 = t81 * t63;
t49 = t81 * t76;
t71 = -t75 * t48 - t62 * t49;
t70 = t76 * MDP(10);
t69 = t75 * MDP(17);
t68 = -t79 - t80;
t36 = -t53 * qJ(5) + t71;
t66 = t62 * t48 - t49 * t75;
t37 = -t51 * qJ(5) - t66;
t67 = t53 * MDP(14) - t51 * MDP(15) + t71 * MDP(17) + t66 * MDP(18) + t36 * MDP(19) - t37 * MDP(20);
t50 = t53 ^ 2;
t1 = [MDP(1) - 0.2e1 * t56 * t70 + t50 * MDP(12) - 0.2e1 * t53 * t51 * MDP(13) + 0.2e1 * (-t36 * t53 - t37 * t51) * MDP(21) + (t36 ^ 2 + t37 ^ 2) * MDP(22) + (t60 ^ 2 + t61 ^ 2) * MDP(4) * pkin(1) ^ 2 + 0.2e1 * t79 * t65 + (0.2e1 * t42 + 0.2e1 * t46 + t73) * t40 + (0.2e1 * t56 * MDP(11) + MDP(5) * t63 + 0.2e1 * MDP(6) * t76) * t63; (-t36 * t51 + t37 * t53) * MDP(22); MDP(4) + (t51 ^ 2 + t50) * MDP(22); t63 * MDP(7) + t76 * MDP(8) + (-t51 * t77 - t58 * t53) * MDP(21) + (t36 * t58 + t37 * t77) * MDP(22) + (-t63 * MDP(10) - MDP(11) * t76) * t55 + t67; t70 - t63 * MDP(11) + (-t51 * t58 + t53 * t77) * MDP(22) + t68; MDP(16) + MDP(9) + (t78 + t72) * t58 + (0.2e1 * t69 + (MDP(22) * t77 - 0.2e1 * MDP(18) - 0.2e1 * MDP(20)) * t62) * pkin(3); (-t53 * MDP(21) + t36 * MDP(22)) * pkin(4) + t67; -t51 * t74 + t68; MDP(16) + (0.2e1 * pkin(4) + t59) * MDP(19) + pkin(4) * t72 + (t69 + (-MDP(18) - MDP(20)) * t62) * pkin(3); MDP(16) + (t78 + t74) * pkin(4); t73 + t80; 0; 0; 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
