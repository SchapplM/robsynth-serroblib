% Calculate joint inertia matrix for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPRP2_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:24
% EndTime: 2019-12-31 17:49:25
% DurationCPUTime: 0.12s
% Computational Cost: add. (189->55), mult. (318->80), div. (0->0), fcn. (300->6), ass. (0->29)
t46 = sin(pkin(7));
t40 = t46 * pkin(1) + qJ(3);
t63 = pkin(6) + t40;
t62 = cos(qJ(4));
t45 = sin(pkin(8));
t49 = sin(qJ(4));
t61 = t49 * t45;
t47 = cos(pkin(8));
t60 = t45 ^ 2 + t47 ^ 2;
t48 = cos(pkin(7));
t42 = -t48 * pkin(1) - pkin(2);
t59 = t42 * MDP(8);
t58 = t45 * MDP(6);
t57 = t47 * MDP(5);
t56 = MDP(15) - MDP(18);
t55 = t62 * t45;
t54 = t60 * MDP(8);
t53 = -pkin(4) * MDP(19) - MDP(16);
t52 = -MDP(14) + t53;
t51 = MDP(19) * qJ(5) - t56;
t38 = -t47 * pkin(3) + t42;
t37 = t49 * t47 + t55;
t36 = -t62 * t47 + t61;
t35 = t37 ^ 2;
t34 = t63 * t47;
t32 = t62 * t34 - t63 * t61;
t31 = t49 * t34 + t63 * t55;
t30 = t36 * pkin(4) - t37 * qJ(5) + t38;
t1 = [MDP(1) + t35 * MDP(9) + (t30 ^ 2 + t31 ^ 2 + t32 ^ 2) * MDP(19) + (t46 ^ 2 + t48 ^ 2) * MDP(4) * pkin(1) ^ 2 + (-0.2e1 * t57 + 0.2e1 * t58 + t59) * t42 + 0.2e1 * (t31 * t37 - t32 * t36) * MDP(17) + (0.2e1 * t60 * MDP(7) + t54 * t40) * t40 + 0.2e1 * (t38 * MDP(15) - t30 * MDP(18)) * t37 + 0.2e1 * (-t37 * MDP(10) + t38 * MDP(14) + t30 * MDP(16)) * t36; (t31 * t36 + t32 * t37) * MDP(19); MDP(4) + t54 + (t36 ^ 2 + t35) * MDP(19); t30 * MDP(19) - t57 + t58 + t59 + t56 * t37 + (MDP(14) + MDP(16)) * t36; 0; MDP(8) + MDP(19); t37 * MDP(11) - t36 * MDP(12) + (-pkin(4) * t37 - t36 * qJ(5)) * MDP(17) + t51 * t32 + t52 * t31; t52 * t36 + t51 * t37; 0; MDP(13) + 0.2e1 * pkin(4) * MDP(16) + 0.2e1 * qJ(5) * MDP(18) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(19); t37 * MDP(17) + t31 * MDP(19); t36 * MDP(19); 0; t53; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
