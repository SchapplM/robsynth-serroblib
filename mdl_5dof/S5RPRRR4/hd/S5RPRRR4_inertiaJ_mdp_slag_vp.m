% Calculate joint inertia matrix for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRRR4_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:22
% EndTime: 2020-01-03 11:52:23
% DurationCPUTime: 0.13s
% Computational Cost: add. (170->53), mult. (287->61), div. (0->0), fcn. (222->8), ass. (0->36)
t56 = sin(qJ(5));
t59 = cos(qJ(5));
t70 = t59 * MDP(16);
t66 = -MDP(17) * t56 + t70;
t63 = 0.2e1 * t66;
t54 = sin(pkin(9));
t79 = pkin(1) * t54;
t78 = pkin(4) * t56;
t60 = cos(qJ(4));
t77 = t60 * pkin(3);
t55 = cos(pkin(9));
t47 = pkin(1) * t55 + pkin(2);
t58 = sin(qJ(3));
t61 = cos(qJ(3));
t43 = t47 * t58 + t61 * t79;
t76 = t60 * t43;
t75 = t56 * MDP(13) + t59 * MDP(14);
t42 = t61 * t47 - t58 * t79;
t41 = pkin(3) + t42;
t57 = sin(qJ(4));
t68 = -t60 * t41 + t43 * t57;
t74 = t68 * MDP(9);
t73 = t42 * MDP(6);
t72 = t43 * MDP(7);
t39 = -t41 * t57 - t76;
t71 = t39 * MDP(10);
t69 = MDP(8) + (MDP(11) * t56 + 0.2e1 * MDP(12) * t59) * t56;
t67 = MDP(5) + t69;
t65 = -MDP(16) * t56 - MDP(17) * t59;
t64 = (-MDP(10) * t57 + MDP(9) * t60) * pkin(3);
t53 = pkin(4) * t59;
t49 = -pkin(4) - t77;
t45 = t49 * t56;
t36 = -pkin(4) + t68;
t35 = t36 * t56;
t1 = [MDP(1) + (t54 ^ 2 + t55 ^ 2) * MDP(4) * pkin(1) ^ 2 - t36 * t63 + 0.2e1 * t71 + 0.2e1 * t73 - 0.2e1 * t72 - 0.2e1 * t74 + t67; 0; MDP(4); t73 - t72 + (-t68 + t77) * MDP(9) + (-t76 + (-pkin(3) - t41) * t57) * MDP(10) + (t45 + t35) * MDP(17) + (-t36 - t49) * t70 + t67; 0; -t49 * t63 + 0.2e1 * t64 + t67; -t74 + t71 + (-t36 * t59 + t53) * MDP(16) + (t35 - t78) * MDP(17) + t69; 0; (-t49 * t59 + t53) * MDP(16) + (t45 - t78) * MDP(17) + t64 + t69; pkin(4) * t63 + t69; t65 * (pkin(8) - t39) + t75; t66; t65 * (pkin(3) * t57 + pkin(8)) + t75; pkin(8) * t65 + t75; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
