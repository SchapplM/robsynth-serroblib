% Calculate joint inertia matrix for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPR6_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:49
% EndTime: 2019-12-31 18:17:50
% DurationCPUTime: 0.08s
% Computational Cost: add. (101->44), mult. (167->49), div. (0->0), fcn. (107->6), ass. (0->24)
t38 = sin(qJ(5));
t40 = cos(qJ(5));
t58 = t38 * MDP(16) + t40 * MDP(17);
t59 = (2 * MDP(9)) + 0.2e1 * t58;
t37 = cos(pkin(8));
t33 = t37 * pkin(1) + pkin(2);
t39 = sin(qJ(3));
t41 = cos(qJ(3));
t36 = sin(pkin(8));
t55 = pkin(1) * t36;
t29 = t41 * t33 - t39 * t55;
t28 = -pkin(3) - t29;
t57 = t28 * MDP(10);
t54 = pkin(3) * MDP(10);
t53 = t29 * MDP(6);
t44 = t39 * t33 + t41 * t55;
t52 = t44 * MDP(7);
t49 = t40 * MDP(16);
t47 = MDP(5) + (MDP(11) * t40 - 0.2e1 * MDP(12) * t38) * t40;
t46 = -t38 * MDP(17) + t49;
t42 = -pkin(3) - pkin(7);
t35 = t40 * MDP(13);
t26 = qJ(4) + t44;
t1 = [MDP(1) + (t36 ^ 2 + t37 ^ 2) * MDP(4) * pkin(1) ^ 2 + 0.2e1 * t53 - 0.2e1 * t52 + t47 + ((2 * MDP(8)) + t57) * t28 + (t26 * MDP(10) + t59) * t26; 0; MDP(4) + MDP(10); t53 - t52 + (-0.2e1 * pkin(3) - t29) * MDP(8) + (0.2e1 * qJ(4) + t44) * MDP(9) + (-t28 * pkin(3) + t26 * qJ(4)) * MDP(10) + t47 + t58 * (qJ(4) + t26); 0; (-(2 * MDP(8)) + t54) * pkin(3) + (MDP(10) * qJ(4) + t59) * qJ(4) + t47; MDP(8) + t57; 0; MDP(8) - t54; MDP(10); -t38 * MDP(14) + t35 + t46 * (-pkin(7) + t28); -t58; t42 * t49 + t35 + (-MDP(17) * t42 - MDP(14)) * t38; t46; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
