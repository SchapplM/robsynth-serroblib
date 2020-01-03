% Calculate joint inertia matrix for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPPRR7_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:54
% EndTime: 2019-12-31 17:59:54
% DurationCPUTime: 0.14s
% Computational Cost: add. (101->56), mult. (188->79), div. (0->0), fcn. (135->6), ass. (0->27)
t44 = sin(qJ(5));
t46 = cos(qJ(5));
t63 = t44 * MDP(20) + t46 * MDP(21);
t43 = cos(pkin(8));
t37 = -t43 * pkin(1) - pkin(2);
t34 = -pkin(6) + t37;
t61 = t34 * t44;
t60 = t34 * t46;
t59 = t44 * t46;
t47 = cos(qJ(4));
t55 = t47 * MDP(14);
t54 = MDP(16) * t59;
t42 = sin(pkin(8));
t35 = t42 * pkin(1) + qJ(3);
t53 = MDP(17) * t46 - MDP(18) * t44;
t52 = t46 * MDP(20) - t44 * MDP(21);
t50 = MDP(13) + t52;
t49 = t44 * MDP(17) + t46 * MDP(18) - pkin(7) * t63;
t45 = sin(qJ(4));
t41 = t47 ^ 2;
t40 = t46 ^ 2;
t39 = t45 ^ 2;
t38 = t44 ^ 2;
t33 = t45 * pkin(4) - t47 * pkin(7) + t35;
t32 = t44 * t33 + t45 * t60;
t31 = t46 * t33 - t45 * t61;
t1 = [MDP(1) + (t35 ^ 2 + t37 ^ 2) * MDP(7) + t39 * MDP(19) + (t42 ^ 2 + t43 ^ 2) * MDP(4) * pkin(1) ^ 2 + (t40 * MDP(15) + MDP(8) - 0.2e1 * t54) * t41 + 0.2e1 * (-MDP(9) + t53) * t47 * t45 + 0.2e1 * t37 * MDP(5) + 0.2e1 * (t31 * t45 - t41 * t61) * MDP(20) + 0.2e1 * (-t32 * t45 - t41 * t60) * MDP(21) + 0.2e1 * (MDP(13) * t45 + MDP(6) + t55) * t35; 0; MDP(4) + MDP(7); t37 * MDP(7) + MDP(5) + t63 * (-t39 - t41); 0; MDP(7); (-t34 * MDP(14) - MDP(11) + t49) * t45 + (MDP(10) + t34 * MDP(13) + MDP(15) * t59 + (-t38 + t40) * MDP(16) + (-pkin(4) * t44 + t60) * MDP(20) + (-pkin(4) * t46 - t61) * MDP(21)) * t47; -t50 * t45 - t55; -t45 * MDP(14) + t50 * t47; t38 * MDP(15) + 0.2e1 * pkin(4) * t52 + MDP(12) + 0.2e1 * t54; t45 * MDP(19) + t31 * MDP(20) - t32 * MDP(21) + t53 * t47; -t63 * t47; -t63 * t45; t49; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
