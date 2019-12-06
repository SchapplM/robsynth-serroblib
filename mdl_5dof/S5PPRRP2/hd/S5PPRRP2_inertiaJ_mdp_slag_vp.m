% Calculate joint inertia matrix for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRRP2_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:09:10
% EndTime: 2019-12-05 15:09:11
% DurationCPUTime: 0.08s
% Computational Cost: add. (83->42), mult. (171->64), div. (0->0), fcn. (143->6), ass. (0->19)
t34 = sin(qJ(4));
t30 = t34 ^ 2;
t36 = cos(qJ(4));
t46 = t36 ^ 2 + t30;
t45 = MDP(12) - MDP(15);
t44 = t46 * MDP(14);
t43 = t46 * MDP(16);
t42 = -pkin(4) * MDP(16) - MDP(13);
t41 = MDP(11) - t42;
t40 = MDP(16) * qJ(5) - t45;
t39 = -t41 * t34 + t40 * t36;
t37 = cos(qJ(3));
t35 = sin(qJ(3));
t33 = cos(pkin(8));
t32 = sin(pkin(8));
t29 = -t36 * pkin(4) - t34 * qJ(5) - pkin(3);
t28 = t37 * t32 + t35 * t33;
t27 = t35 * t32 - t37 * t33;
t1 = [MDP(1) + (t32 ^ 2 + t33 ^ 2) * MDP(2) + (t46 * t28 ^ 2 + t27 ^ 2) * MDP(16); 0; MDP(2) + t43; (pkin(6) * t43 - MDP(5) + t44) * t28 + (t29 * MDP(16) - MDP(4) + (-MDP(11) - MDP(13)) * t36 + t45 * t34) * t27; 0; MDP(3) + t30 * MDP(6) + (t46 * pkin(6) ^ 2 + t29 ^ 2) * MDP(16) + 0.2e1 * pkin(6) * t44 + 0.2e1 * (pkin(3) * MDP(11) - t29 * MDP(13)) * t36 + 0.2e1 * (-pkin(3) * MDP(12) - t29 * MDP(15) + t36 * MDP(7)) * t34; t39 * t28; t40 * t34 + t41 * t36; t34 * MDP(8) + t36 * MDP(9) + (-t34 * pkin(4) + t36 * qJ(5)) * MDP(14) + t39 * pkin(6); MDP(10) + 0.2e1 * pkin(4) * MDP(13) + 0.2e1 * qJ(5) * MDP(15) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(16); t34 * t28 * MDP(16); -t36 * MDP(16); (pkin(6) * MDP(16) + MDP(14)) * t34; t42; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
