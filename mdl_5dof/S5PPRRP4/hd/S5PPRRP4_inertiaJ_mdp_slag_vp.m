% Calculate joint inertia matrix for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRRP4_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:56:29
% EndTime: 2021-01-15 14:56:30
% DurationCPUTime: 0.09s
% Computational Cost: add. (64->43), mult. (125->55), div. (0->0), fcn. (87->4), ass. (0->20)
t40 = qJ(5) + pkin(6);
t27 = sin(qJ(4));
t24 = t27 ^ 2;
t29 = cos(qJ(4));
t39 = t29 ^ 2 + t24;
t38 = MDP(16) * pkin(4);
t22 = -t29 * pkin(4) - pkin(3);
t37 = t22 * MDP(16);
t36 = t27 * MDP(12);
t23 = t27 * MDP(14);
t35 = t29 * MDP(13);
t34 = -MDP(12) - MDP(14);
t33 = MDP(13) + t38;
t20 = t40 * t27;
t21 = t40 * t29;
t32 = t20 * t27 + t21 * t29;
t31 = -MDP(11) - t33;
t30 = cos(qJ(3));
t28 = sin(qJ(3));
t1 = [t39 * MDP(16) + MDP(1) + MDP(2); 0; MDP(2) + (t39 * t28 ^ 2 + t30 ^ 2) * MDP(16); (t29 * t20 - t27 * t21) * MDP(16); (t39 * MDP(15) + t32 * MDP(16) - MDP(5)) * t28 + (-t37 + MDP(4) + (MDP(11) + MDP(13)) * t29 + t34 * t27) * t30; MDP(3) + t24 * MDP(6) + 0.2e1 * t27 * t29 * MDP(7) + 0.2e1 * t32 * MDP(15) + (t20 ^ 2 + t21 ^ 2) * MDP(16) + (0.2e1 * t23 - 0.2e1 * t35 + t37) * t22 + 0.2e1 * (t29 * MDP(11) - t36) * pkin(3); t31 * t29 + t23 + t36; (t31 * t27 + t34 * t29) * t28; -t21 * MDP(14) + (-MDP(12) * pkin(6) + MDP(9)) * t29 - t33 * t20 + (-MDP(11) * pkin(6) - MDP(15) * pkin(4) + MDP(8)) * t27; MDP(10) + (0.2e1 * MDP(13) + t38) * pkin(4); 0; -t30 * MDP(16); t23 - t35 + t37; 0; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
