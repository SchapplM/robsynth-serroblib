% Calculate joint inertia matrix for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RRPP3_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:51
% EndTime: 2019-12-31 16:57:51
% DurationCPUTime: 0.10s
% Computational Cost: add. (130->49), mult. (237->74), div. (0->0), fcn. (204->4), ass. (0->21)
t44 = cos(qJ(2));
t53 = 0.2e1 * t44;
t52 = -qJ(3) - pkin(5);
t41 = sin(pkin(6));
t42 = cos(pkin(6));
t43 = sin(qJ(2));
t32 = t41 * t43 - t42 * t44;
t33 = t41 * t44 + t42 * t43;
t40 = -t44 * pkin(2) - pkin(1);
t27 = t32 * pkin(3) - t33 * qJ(4) + t40;
t51 = t27 * MDP(16);
t50 = t32 * MDP(13);
t49 = t33 * MDP(15);
t35 = t52 * t44;
t47 = t52 * t43;
t29 = -t41 * t35 - t42 * t47;
t31 = -t42 * t35 + t41 * t47;
t48 = t29 ^ 2 + t31 ^ 2;
t38 = t42 * pkin(2) + pkin(3);
t36 = t41 * pkin(2) + qJ(4);
t1 = [MDP(1) + pkin(1) * MDP(9) * t53 + (t40 ^ 2 + t48) * MDP(12) + t48 * MDP(16) + (-0.2e1 * t49 + 0.2e1 * t50 + t51) * t27 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t43 + MDP(5) * t53) * t43 + 0.2e1 * (MDP(11) + MDP(14)) * (t29 * t33 - t31 * t32); t43 * MDP(6) + t44 * MDP(7) - t29 * MDP(13) + (-t36 * t32 - t38 * t33) * MDP(14) + t31 * MDP(15) + (-t29 * t38 + t31 * t36) * MDP(16) + (-t44 * MDP(10) - t43 * MDP(9)) * pkin(5) + ((-t32 * t41 - t33 * t42) * MDP(11) + (-t29 * t42 + t31 * t41) * MDP(12)) * pkin(2); MDP(8) + (t36 ^ 2 + t38 ^ 2) * MDP(16) + (t41 ^ 2 + t42 ^ 2) * MDP(12) * pkin(2) ^ 2 + 0.2e1 * t38 * MDP(13) + 0.2e1 * t36 * MDP(15); t40 * MDP(12) - t49 + t50 + t51; 0; MDP(12) + MDP(16); t33 * MDP(14) + t29 * MDP(16); -t38 * MDP(16) - MDP(13); 0; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
