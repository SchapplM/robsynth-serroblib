% Calculate joint inertia matrix for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR10_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:56
% EndTime: 2019-12-31 17:11:56
% DurationCPUTime: 0.13s
% Computational Cost: add. (94->59), mult. (188->79), div. (0->0), fcn. (124->4), ass. (0->32)
t47 = cos(qJ(2));
t65 = 0.2e1 * t47;
t48 = -pkin(2) - pkin(6);
t64 = pkin(3) + pkin(5);
t39 = t64 * t47;
t63 = t39 * t47;
t44 = sin(qJ(4));
t46 = cos(qJ(4));
t62 = t44 * t46;
t61 = pkin(2) * MDP(14);
t60 = t44 * MDP(20);
t59 = t46 * MDP(21);
t58 = pkin(5) ^ 2 * MDP(14);
t57 = MDP(14) * qJ(3);
t56 = MDP(16) * t62;
t45 = sin(qJ(2));
t55 = -t45 * qJ(3) - pkin(1);
t54 = MDP(12) - t61;
t53 = -MDP(17) * t44 - MDP(18) * t46;
t52 = t46 * MDP(20) - t44 * MDP(21);
t51 = MDP(11) + t52;
t50 = (MDP(20) * t48 + MDP(17)) * t46 + (-MDP(21) * t48 - MDP(18)) * t44;
t43 = t47 ^ 2;
t42 = t46 ^ 2;
t41 = t45 ^ 2;
t40 = t44 ^ 2;
t38 = t64 * t45;
t37 = -t47 * pkin(2) + t55;
t36 = t48 * t47 + t55;
t35 = t46 * t36 + t44 * t38;
t34 = -t44 * t36 + t46 * t38;
t1 = [pkin(1) * MDP(9) * t65 + MDP(1) + (MDP(12) * t65 + MDP(14) * t37) * t37 + (t40 * MDP(15) + 0.2e1 * t56 + t58) * t43 + (MDP(19) + MDP(4) + t58) * t41 + 0.2e1 * (-pkin(1) * MDP(10) - t37 * MDP(13) + (MDP(5) + t53) * t47) * t45 + 0.2e1 * (t34 * t45 + t46 * t63) * MDP(20) + 0.2e1 * (-t35 * t45 - t44 * t63) * MDP(21) + 0.2e1 * (t41 + t43) * MDP(11) * pkin(5); (t59 + t60) * t39 + (-pkin(2) * MDP(11) + MDP(6) + t50) * t45 + (MDP(7) - MDP(15) * t62 + (t40 - t42) * MDP(16) + t51 * qJ(3)) * t47 + ((-MDP(10) + MDP(13) + t57) * t47 + (-MDP(9) + t54) * t45) * pkin(5); -0.2e1 * t56 + t42 * MDP(15) + MDP(8) + (-0.2e1 * MDP(12) + t61) * pkin(2) + (0.2e1 * MDP(13) + t57 + 0.2e1 * t59 + 0.2e1 * t60) * qJ(3); (MDP(14) * pkin(5) + t51) * t45; t54; MDP(14); t45 * MDP(19) + t34 * MDP(20) - t35 * MDP(21) + t53 * t47; t50; t52; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
