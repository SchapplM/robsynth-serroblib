% Calculate joint inertia matrix for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPPR3_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:02
% EndTime: 2019-12-31 17:44:02
% DurationCPUTime: 0.09s
% Computational Cost: add. (93->42), mult. (157->59), div. (0->0), fcn. (120->6), ass. (0->25)
t54 = MDP(8) + MDP(12);
t46 = sin(pkin(8));
t48 = cos(pkin(8));
t57 = t46 ^ 2 + t48 ^ 2;
t63 = t54 * t57;
t49 = cos(pkin(7));
t43 = -t49 * pkin(1) - pkin(2);
t53 = t46 * qJ(4) - t43;
t62 = 0.2e1 * (pkin(3) + pkin(4)) * t48 + 0.2e1 * t53;
t61 = -0.2e1 * t48;
t47 = sin(pkin(7));
t42 = t47 * pkin(1) + qJ(3);
t60 = -pkin(6) + t42;
t50 = sin(qJ(5));
t51 = cos(qJ(5));
t36 = t46 * t50 + t48 * t51;
t34 = t36 * MDP(18);
t37 = t46 * t51 - t48 * t50;
t59 = -t37 * MDP(19) - t34;
t56 = t43 * MDP(8);
t31 = -t48 * pkin(3) - t53;
t55 = t31 * MDP(12);
t33 = t60 * t48;
t32 = t60 * t46;
t1 = [MDP(1) + t34 * t62 + (t47 ^ 2 + t49 ^ 2) * MDP(4) * pkin(1) ^ 2 + (MDP(5) * t61 + 0.2e1 * t46 * MDP(6) + t56) * t43 + (-0.2e1 * t46 * MDP(11) + MDP(9) * t61 + t55) * t31 + (MDP(13) * t37 - 0.2e1 * t36 * MDP(14) + MDP(19) * t62) * t37 + (0.2e1 * (MDP(10) + MDP(7)) * t57 + t63 * t42) * t42; 0; MDP(4) + t63; t55 + t56 + (-MDP(5) - MDP(9)) * t48 + (-MDP(11) + MDP(6)) * t46 + t59; 0; t54; (MDP(12) * t42 + MDP(10)) * t46; -t48 * MDP(12); 0; MDP(12); t37 * MDP(15) - t36 * MDP(16) + (t51 * t32 - t50 * t33) * MDP(18) + (-t50 * t32 - t51 * t33) * MDP(19); t59; 0; t51 * MDP(18) - t50 * MDP(19); MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
