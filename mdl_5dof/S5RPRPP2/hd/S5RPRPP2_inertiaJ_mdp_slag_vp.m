% Calculate joint inertia matrix for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPRPP2_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:12
% EndTime: 2019-12-31 18:11:13
% DurationCPUTime: 0.13s
% Computational Cost: add. (133->66), mult. (199->83), div. (0->0), fcn. (123->4), ass. (0->24)
t46 = sin(qJ(3));
t47 = cos(qJ(3));
t60 = t47 * MDP(16) + t46 * MDP(17);
t40 = t46 * qJ(4);
t59 = t47 * pkin(3) + t40;
t42 = t46 ^ 2;
t58 = t47 ^ 2 + t42;
t57 = pkin(3) * MDP(15);
t56 = t47 * qJ(4);
t44 = sin(pkin(7));
t36 = t44 * pkin(1) + pkin(6);
t55 = -qJ(5) + t36;
t54 = -MDP(11) + MDP(14);
t53 = MDP(15) + MDP(19);
t45 = cos(pkin(7));
t37 = -t45 * pkin(1) - pkin(2);
t48 = pkin(3) + pkin(4);
t52 = t48 * MDP(19) + MDP(12);
t31 = t37 - t59;
t50 = qJ(4) ^ 2;
t33 = t55 * t47;
t32 = t55 * t46;
t30 = t47 * pkin(4) - t31;
t1 = [MDP(1) + t42 * MDP(5) + (t58 * t36 ^ 2 + t31 ^ 2) * MDP(15) + (t30 ^ 2 + t32 ^ 2 + t33 ^ 2) * MDP(19) + (t44 ^ 2 + t45 ^ 2) * MDP(4) * pkin(1) ^ 2 + 0.2e1 * t58 * MDP(13) * t36 + 0.2e1 * (-t37 * MDP(10) - t31 * MDP(12) + t30 * MDP(16) - MDP(18) * t33) * t47 + 0.2e1 * (t37 * MDP(11) - t31 * MDP(14) + t30 * MDP(17) - MDP(18) * t32 + t47 * MDP(6)) * t46; (-t32 * t47 + t33 * t46) * MDP(19); t53 * t58 + MDP(4); t46 * MDP(7) + t47 * MDP(8) + (-t46 * pkin(3) + t56) * MDP(13) - t32 * MDP(16) + t33 * MDP(17) + (t48 * t46 - t56) * MDP(18) + (t33 * qJ(4) - t32 * t48) * MDP(19) + ((qJ(4) * MDP(15) + t54) * t47 + (-MDP(10) - MDP(12) - t57) * t46) * t36; t59 * MDP(15) + t40 * MDP(19) + (MDP(10) + t52) * t47 + t54 * t46 + t60; MDP(9) + 0.2e1 * pkin(3) * MDP(12) + (pkin(3) ^ 2 + t50) * MDP(15) + 0.2e1 * t48 * MDP(16) + (t48 ^ 2 + t50) * MDP(19) + 0.2e1 * (MDP(14) + MDP(17)) * qJ(4); t32 * MDP(19) + (MDP(15) * t36 + MDP(13) - MDP(18)) * t46; -t53 * t47; -MDP(16) - t52 - t57; t53; t30 * MDP(19) + t60; 0; 0; 0; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
