% Calculate joint inertia matrix for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RPPR7_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:41
% EndTime: 2019-12-31 16:41:41
% DurationCPUTime: 0.05s
% Computational Cost: add. (59->32), mult. (98->49), div. (0->0), fcn. (72->4), ass. (0->17)
t30 = sin(pkin(6));
t39 = 0.2e1 * t30 * pkin(3) + (2 * qJ(2));
t32 = -pkin(1) - qJ(3);
t38 = -pkin(5) + t32;
t31 = cos(pkin(6));
t25 = t30 ^ 2 + t31 ^ 2;
t33 = sin(qJ(4));
t34 = cos(qJ(4));
t21 = t34 * t30 + t33 * t31;
t37 = t21 * MDP(16);
t36 = t30 * MDP(7) + t31 * MDP(8);
t35 = qJ(2) ^ 2;
t24 = t38 * t31;
t23 = t38 * t30;
t22 = -t33 * t30 + t34 * t31;
t20 = t25 * t32;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2 + t35) * MDP(6)) - 0.2e1 * t20 * MDP(9) + (t25 * t32 ^ 2 + t35) * MDP(10) + t37 * t39 + (MDP(11) * t22 - 0.2e1 * t21 * MDP(12) + MDP(17) * t39) * t22 + 0.2e1 * (MDP(5) + t36) * qJ(2); t20 * MDP(10) - (pkin(1) * MDP(6)) - t25 * MDP(9) + MDP(4); t25 * MDP(10) + MDP(6); qJ(2) * MDP(10) + t22 * MDP(17) + t36 + t37; 0; MDP(10); t22 * MDP(13) - t21 * MDP(14) + (-t33 * t23 + t34 * t24) * MDP(16) + (-t34 * t23 - t33 * t24) * MDP(17); t22 * MDP(16) - t21 * MDP(17); 0; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
