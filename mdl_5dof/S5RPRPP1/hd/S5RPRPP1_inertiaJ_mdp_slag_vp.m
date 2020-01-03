% Calculate joint inertia matrix for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPP1_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:14
% EndTime: 2019-12-31 18:09:15
% DurationCPUTime: 0.16s
% Computational Cost: add. (235->62), mult. (395->96), div. (0->0), fcn. (376->6), ass. (0->27)
t79 = MDP(13) + MDP(17);
t58 = sin(pkin(8));
t60 = cos(pkin(8));
t62 = sin(qJ(3));
t63 = cos(qJ(3));
t49 = t58 * t62 - t60 * t63;
t51 = t58 * t63 + t60 * t62;
t61 = cos(pkin(7));
t57 = -t61 * pkin(1) - pkin(2);
t52 = -t63 * pkin(3) + t57;
t40 = t49 * pkin(4) - t51 * qJ(5) + t52;
t78 = t40 * MDP(17);
t77 = t49 * MDP(14);
t76 = t51 * MDP(16);
t75 = t63 * MDP(10);
t59 = sin(pkin(7));
t72 = t59 * pkin(1) + pkin(6);
t68 = qJ(4) + t72;
t47 = t68 * t63;
t66 = t68 * t62;
t42 = t58 * t47 + t60 * t66;
t44 = t60 * t47 - t58 * t66;
t74 = t42 ^ 2 + t44 ^ 2;
t67 = t76 - t77;
t55 = t60 * pkin(3) + pkin(4);
t53 = t58 * pkin(3) + qJ(5);
t1 = [MDP(1) - 0.2e1 * t57 * t75 + (t52 ^ 2 + t74) * MDP(13) + t74 * MDP(17) + (t59 ^ 2 + t61 ^ 2) * MDP(4) * pkin(1) ^ 2 + (-0.2e1 * t76 + 0.2e1 * t77 + t78) * t40 + (0.2e1 * t57 * MDP(11) + MDP(5) * t62 + 0.2e1 * t63 * MDP(6)) * t62 + 0.2e1 * (MDP(12) + MDP(15)) * (t42 * t51 - t44 * t49); t79 * (t42 * t49 + t44 * t51); MDP(4) + t79 * (t49 ^ 2 + t51 ^ 2); -t42 * MDP(14) + (-t53 * t49 - t55 * t51) * MDP(15) + t44 * MDP(16) + (-t42 * t55 + t44 * t53) * MDP(17) + (-t72 * MDP(11) + MDP(8)) * t63 + (-t72 * MDP(10) + MDP(7)) * t62 + ((-t49 * t58 - t51 * t60) * MDP(12) + (-t42 * t60 + t44 * t58) * MDP(13)) * pkin(3); t75 - t62 * MDP(11) + (-t49 * t55 + t51 * t53) * MDP(17) + (-t49 * t60 + t51 * t58) * MDP(13) * pkin(3) + t67; MDP(9) + (t53 ^ 2 + t55 ^ 2) * MDP(17) + (t58 ^ 2 + t60 ^ 2) * MDP(13) * pkin(3) ^ 2 + 0.2e1 * t55 * MDP(14) + 0.2e1 * t53 * MDP(16); t52 * MDP(13) - t67 + t78; 0; 0; t79; t51 * MDP(15) + t42 * MDP(17); t49 * MDP(17); -t55 * MDP(17) - MDP(14); 0; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
