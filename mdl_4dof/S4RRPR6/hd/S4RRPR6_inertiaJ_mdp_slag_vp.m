% Calculate joint inertia matrix for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S4RRPR6_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:55
% EndTime: 2019-12-31 17:04:55
% DurationCPUTime: 0.08s
% Computational Cost: add. (151->49), mult. (289->80), div. (0->0), fcn. (293->6), ass. (0->27)
t56 = sin(pkin(7));
t57 = cos(pkin(7));
t59 = sin(qJ(2));
t61 = cos(qJ(2));
t47 = -t56 * t59 + t57 * t61;
t55 = -t61 * pkin(2) - pkin(1);
t70 = -0.2e1 * t47 * pkin(3) + 0.2e1 * t55;
t69 = 0.2e1 * t61;
t68 = pkin(2) * t56;
t67 = -qJ(3) - pkin(5);
t52 = t67 * t59;
t53 = t67 * t61;
t42 = t56 * t52 - t57 * t53;
t48 = t56 * t61 + t57 * t59;
t58 = sin(qJ(4));
t60 = cos(qJ(4));
t39 = -t60 * t47 + t58 * t48;
t66 = t39 * MDP(18);
t54 = t57 * pkin(2) + pkin(3);
t65 = (t60 * t54 - t58 * t68) * MDP(18);
t64 = (-t58 * t54 - t60 * t68) * MDP(19);
t41 = t57 * t52 + t56 * t53;
t35 = -t48 * pkin(6) + t41;
t36 = t47 * pkin(6) + t42;
t40 = t58 * t47 + t60 * t48;
t63 = t40 * MDP(15) - t39 * MDP(16) + (t60 * t35 - t58 * t36) * MDP(18) + (-t58 * t35 - t60 * t36) * MDP(19);
t1 = [MDP(1) + pkin(1) * MDP(9) * t69 + 0.2e1 * (-t41 * t48 + t42 * t47) * MDP(11) + (t41 ^ 2 + t42 ^ 2 + t55 ^ 2) * MDP(12) + t66 * t70 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t59 + MDP(5) * t69) * t59 + (MDP(13) * t40 - 0.2e1 * t39 * MDP(14) + MDP(19) * t70) * t40; t59 * MDP(6) + t61 * MDP(7) + (-t61 * MDP(10) - t59 * MDP(9)) * pkin(5) + ((t47 * t56 - t48 * t57) * MDP(11) + (t41 * t57 + t42 * t56) * MDP(12)) * pkin(2) + t63; MDP(8) + MDP(17) + (t56 ^ 2 + t57 ^ 2) * MDP(12) * pkin(2) ^ 2 + 0.2e1 * t65 + 0.2e1 * t64; t55 * MDP(12) + t40 * MDP(19) + t66; 0; MDP(12); t63; MDP(17) + t64 + t65; 0; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
