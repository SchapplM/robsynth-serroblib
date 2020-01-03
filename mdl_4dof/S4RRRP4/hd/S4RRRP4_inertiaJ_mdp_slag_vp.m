% Calculate joint inertia matrix for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S4RRRP4_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:42
% EndTime: 2019-12-31 17:15:42
% DurationCPUTime: 0.08s
% Computational Cost: add. (111->45), mult. (206->74), div. (0->0), fcn. (182->4), ass. (0->21)
t41 = cos(qJ(2));
t49 = 0.2e1 * t41;
t48 = pkin(5) + pkin(6);
t38 = sin(qJ(3));
t47 = pkin(2) * t38;
t40 = cos(qJ(3));
t46 = t40 * MDP(16);
t45 = -t41 * pkin(2) - pkin(1);
t39 = sin(qJ(2));
t33 = t48 * t39;
t34 = t48 * t41;
t44 = -t40 * t33 - t38 * t34;
t30 = t38 * t39 - t40 * t41;
t31 = t38 * t41 + t40 * t39;
t42 = t38 * t33 - t40 * t34;
t43 = t31 * MDP(13) - t30 * MDP(14) + t44 * MDP(16) + t42 * MDP(17);
t36 = t40 * pkin(2) + pkin(3);
t27 = t30 * pkin(3) + t45;
t24 = -t30 * qJ(4) - t42;
t23 = -t31 * qJ(4) + t44;
t1 = [MDP(1) + pkin(1) * MDP(9) * t49 + t31 ^ 2 * MDP(11) - 0.2e1 * t31 * t30 * MDP(12) + 0.2e1 * (-t23 * t31 - t24 * t30) * MDP(18) + (t23 ^ 2 + t24 ^ 2 + t27 ^ 2) * MDP(19) + 0.2e1 * (t30 * MDP(16) + t31 * MDP(17)) * t45 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t39 + MDP(5) * t49) * t39; t39 * MDP(6) + t41 * MDP(7) + (-t30 * t47 - t36 * t31) * MDP(18) + (t23 * t36 + t24 * t47) * MDP(19) + (-t41 * MDP(10) - t39 * MDP(9)) * pkin(5) + t43; t36 ^ 2 * MDP(19) + MDP(15) + MDP(8) + (0.2e1 * t46 + (MDP(19) * t47 - 0.2e1 * MDP(17)) * t38) * pkin(2); (-t31 * MDP(18) + t23 * MDP(19)) * pkin(3) + t43; t36 * pkin(3) * MDP(19) + MDP(15) + (-t38 * MDP(17) + t46) * pkin(2); pkin(3) ^ 2 * MDP(19) + MDP(15); t27 * MDP(19); 0; 0; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
