% Calculate joint inertia matrix for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RRPP4_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:19
% EndTime: 2019-12-31 16:59:19
% DurationCPUTime: 0.09s
% Computational Cost: add. (82->53), mult. (134->67), div. (0->0), fcn. (68->2), ass. (0->15)
t30 = pkin(2) + pkin(3);
t38 = pkin(5) - qJ(4);
t28 = sin(qJ(2));
t26 = t28 ^ 2;
t29 = cos(qJ(2));
t37 = t29 ^ 2 + t26;
t36 = t29 * qJ(3);
t35 = t28 * qJ(3) + pkin(1);
t34 = -pkin(2) * MDP(14) - MDP(11);
t32 = qJ(3) ^ 2;
t25 = t38 * t29;
t24 = t38 * t28;
t23 = -t29 * pkin(2) - t35;
t22 = t30 * t29 + t35;
t1 = [MDP(1) + t26 * MDP(4) + (t37 * pkin(5) ^ 2 + t23 ^ 2) * MDP(14) + (t22 ^ 2 + t24 ^ 2 + t25 ^ 2) * MDP(18) + 0.2e1 * t37 * MDP(12) * pkin(5) + 0.2e1 * (-t23 * MDP(11) + t22 * MDP(15) - MDP(17) * t25 + pkin(1) * MDP(9)) * t29 + 0.2e1 * (-pkin(1) * MDP(10) - t23 * MDP(13) + t22 * MDP(16) - MDP(17) * t24 + t29 * MDP(5)) * t28; t28 * MDP(6) + t29 * MDP(7) + (-t28 * pkin(2) + t36) * MDP(12) - t24 * MDP(15) + t25 * MDP(16) + (t30 * t28 - t36) * MDP(17) + (t25 * qJ(3) - t24 * t30) * MDP(18) + ((qJ(3) * MDP(14) - MDP(10) + MDP(13)) * t29 + (-MDP(9) + t34) * t28) * pkin(5); MDP(8) + 0.2e1 * pkin(2) * MDP(11) + (pkin(2) ^ 2 + t32) * MDP(14) + 0.2e1 * t30 * MDP(15) + (t30 ^ 2 + t32) * MDP(18) + 0.2e1 * (MDP(13) + MDP(16)) * qJ(3); t24 * MDP(18) + (MDP(14) * pkin(5) + MDP(12) - MDP(17)) * t28; -t30 * MDP(18) - MDP(15) + t34; MDP(14) + MDP(18); t29 * MDP(15) + t28 * MDP(16) + t22 * MDP(18); 0; 0; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
