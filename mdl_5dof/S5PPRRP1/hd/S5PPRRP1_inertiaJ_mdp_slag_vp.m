% Calculate joint inertia matrix for
% S5PPRRP1
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
%   see S5PPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRRP1_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:48:17
% EndTime: 2021-01-15 14:48:18
% DurationCPUTime: 0.11s
% Computational Cost: add. (84->44), mult. (168->62), div. (0->0), fcn. (149->6), ass. (0->24)
t34 = sin(qJ(4));
t41 = MDP(12) + MDP(14);
t48 = t41 * t34;
t47 = -qJ(5) - pkin(6);
t30 = t34 ^ 2;
t36 = cos(qJ(4));
t46 = t36 ^ 2 + t30;
t45 = MDP(16) * pkin(4);
t29 = -t36 * pkin(4) - pkin(3);
t44 = t29 * MDP(16);
t43 = t34 * MDP(14);
t42 = t36 * MDP(13);
t40 = MDP(13) + t45;
t27 = t47 * t34;
t28 = t47 * t36;
t39 = -t27 * t34 - t28 * t36;
t38 = MDP(11) + t40;
t37 = cos(qJ(3));
t35 = sin(qJ(3));
t33 = cos(pkin(8));
t32 = sin(pkin(8));
t26 = t37 * t32 + t35 * t33;
t25 = t35 * t32 - t37 * t33;
t1 = [MDP(1) + (t32 ^ 2 + t33 ^ 2) * MDP(2) + (t46 * t26 ^ 2 + t25 ^ 2) * MDP(16); 0; t46 * MDP(16) + MDP(2); (t46 * MDP(15) + t39 * MDP(16) - MDP(5)) * t26 + (t44 - MDP(4) + (-MDP(11) - MDP(13)) * t36 + t48) * t25; (t36 * t27 - t34 * t28) * MDP(16); MDP(3) + t30 * MDP(6) + 0.2e1 * t34 * t36 * MDP(7) + 0.2e1 * t39 * MDP(15) + (t27 ^ 2 + t28 ^ 2) * MDP(16) + (-0.2e1 * t42 + 0.2e1 * t43 + t44) * t29 + 0.2e1 * (t36 * MDP(11) - t34 * MDP(12)) * pkin(3); (-t38 * t34 - t41 * t36) * t26; t38 * t36 - t48; t28 * MDP(14) + (-MDP(12) * pkin(6) + MDP(9)) * t36 + t40 * t27 + (-MDP(11) * pkin(6) - MDP(15) * pkin(4) + MDP(8)) * t34; MDP(10) + (0.2e1 * MDP(13) + t45) * pkin(4); t25 * MDP(16); 0; -t42 + t43 + t44; 0; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
