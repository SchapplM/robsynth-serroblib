% Calculate joint inertia matrix for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP12_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:21
% EndTime: 2019-12-31 18:57:22
% DurationCPUTime: 0.21s
% Computational Cost: add. (192->94), mult. (342->130), div. (0->0), fcn. (260->4), ass. (0->41)
t56 = sin(qJ(3));
t80 = 0.2e1 * t56;
t79 = 2 * MDP(21);
t78 = (pkin(1) * MDP(6));
t55 = sin(qJ(4));
t57 = cos(qJ(4));
t77 = t55 * t57;
t59 = -pkin(1) - pkin(6);
t76 = t55 * t59;
t75 = t57 * t59;
t74 = -qJ(5) - pkin(7);
t51 = t55 ^ 2;
t53 = t57 ^ 2;
t73 = t51 + t53;
t52 = t56 ^ 2;
t58 = cos(qJ(3));
t54 = t58 ^ 2;
t72 = -t52 - t54;
t71 = MDP(22) * pkin(4);
t70 = qJ(5) * t58;
t47 = t56 * pkin(3) - t58 * pkin(7) + qJ(2);
t45 = t57 * t47;
t41 = -t57 * t70 + t45 + (pkin(4) - t76) * t56;
t69 = t41 * MDP(22);
t50 = -t57 * pkin(4) - pkin(3);
t68 = t50 * MDP(22);
t67 = t55 * MDP(17);
t66 = t57 * MDP(20);
t65 = t58 * MDP(22);
t64 = t56 * t75;
t63 = MDP(15) * t77;
t62 = -MDP(21) * pkin(4) + MDP(16);
t48 = t74 * t55;
t49 = t74 * t57;
t61 = -t48 * t55 - t49 * t57;
t60 = t57 * MDP(19) - t55 * MDP(20);
t46 = (pkin(4) * t55 - t59) * t58;
t44 = t55 * t47 + t64;
t43 = -t56 * t76 + t45;
t42 = t64 + (t47 - t70) * t55;
t1 = [MDP(1) + t52 * MDP(18) + (t41 ^ 2 + t42 ^ 2 + t46 ^ 2) * MDP(22) + ((-2 * MDP(4) + t78) * pkin(1)) + (MDP(12) * t80 + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + (t53 * MDP(14) + MDP(7) - 0.2e1 * t63) * t54 + 0.2e1 * (t43 * t56 - t54 * t76) * MDP(19) + 0.2e1 * (-t44 * t56 - t54 * t75) * MDP(20) + (0.2e1 * qJ(2) * MDP(13) + (-t41 * t57 - t42 * t55) * t79 + (t57 * MDP(16) - MDP(8) - t67) * t80) * t58; -t46 * t65 - t78 + MDP(4) + (t42 * t56 * MDP(22) + t72 * MDP(20)) * t57 + (t72 * MDP(19) - t56 * t69) * t55; MDP(6) + (t73 * t52 + t54) * MDP(22); (-t41 * t55 + t42 * t57) * MDP(21) + (t41 * t48 - t42 * t49 + t46 * t50) * MDP(22) + (-t59 * MDP(13) + t55 * MDP(16) + t57 * MDP(17) - MDP(10) + (-t55 * MDP(19) - t66) * pkin(7)) * t56 + (MDP(9) + t59 * MDP(12) + MDP(14) * t77 + (-t51 + t53) * MDP(15) + (-pkin(3) * t55 + t75) * MDP(19) + (-pkin(3) * t57 - t76) * MDP(20) + (-t48 * t57 + t49 * t55) * MDP(21)) * t58; (MDP(12) + t60 - t68) * t58 + (t73 * MDP(21) + t61 * MDP(22) - MDP(13)) * t56; MDP(11) + t51 * MDP(14) + 0.2e1 * t63 + t61 * t79 + (t48 ^ 2 + t49 ^ 2 + t50 ^ 2) * MDP(22) + 0.2e1 * t60 * pkin(3); pkin(4) * t69 + t56 * MDP(18) + t43 * MDP(19) - t44 * MDP(20) + (t62 * t57 - t67) * t58; (-t66 + (-MDP(19) - t71) * t55) * t56; t48 * t71 + (-MDP(20) * pkin(7) + MDP(17)) * t57 + (-MDP(19) * pkin(7) + t62) * t55; MDP(22) * pkin(4) ^ 2 + MDP(18); t46 * MDP(22); -t65; t68; 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
