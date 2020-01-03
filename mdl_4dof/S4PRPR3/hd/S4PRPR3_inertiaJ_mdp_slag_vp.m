% Calculate joint inertia matrix for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRPR3_inertiaJ_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:55
% EndTime: 2019-12-31 16:20:55
% DurationCPUTime: 0.08s
% Computational Cost: add. (40->24), mult. (82->38), div. (0->0), fcn. (68->4), ass. (0->18)
t26 = cos(pkin(7));
t38 = -0.2e1 * t26 * pkin(3) - (2 * pkin(2));
t37 = pkin(2) * MDP(8);
t36 = pkin(5) + qJ(3);
t25 = sin(pkin(7));
t35 = t25 ^ 2 + t26 ^ 2;
t34 = t25 * MDP(6);
t33 = t26 * MDP(5);
t27 = sin(qJ(4));
t28 = cos(qJ(4));
t17 = t27 * t25 - t28 * t26;
t32 = t17 * MDP(14);
t31 = t35 * MDP(8);
t18 = t28 * t25 + t27 * t26;
t30 = -t18 * MDP(15) - t32;
t20 = t36 * t26;
t19 = t36 * t25;
t1 = [MDP(1) + t31; 0; t32 * t38 + MDP(2) + (0.2e1 * t33 - 0.2e1 * t34 + t37) * pkin(2) + (-0.2e1 * t17 * MDP(10) + MDP(15) * t38 + MDP(9) * t18) * t18 + (0.2e1 * t35 * MDP(7) + t31 * qJ(3)) * qJ(3); 0; -t30 - t33 + t34 - t37; MDP(8); t30; t18 * MDP(11) - t17 * MDP(12) + (-t28 * t19 - t27 * t20) * MDP(14) + (t27 * t19 - t28 * t20) * MDP(15); 0; MDP(13);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
