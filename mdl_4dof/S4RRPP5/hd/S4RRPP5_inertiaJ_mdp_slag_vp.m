% Calculate joint inertia matrix for
% S4RRPP5
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
%   see S4RRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RRPP5_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:37
% EndTime: 2019-12-31 17:00:38
% DurationCPUTime: 0.08s
% Computational Cost: add. (80->53), mult. (133->67), div. (0->0), fcn. (66->2), ass. (0->15)
t39 = pkin(3) + pkin(5);
t30 = pkin(2) + qJ(4);
t32 = sin(qJ(2));
t28 = t32 ^ 2;
t33 = cos(qJ(2));
t38 = t33 ^ 2 + t28;
t37 = -t32 * qJ(3) - pkin(1);
t36 = -pkin(2) * MDP(14) + MDP(12);
t34 = qJ(3) ^ 2;
t27 = t33 * qJ(3);
t26 = t39 * t33;
t25 = t39 * t32;
t24 = -t33 * pkin(2) + t37;
t23 = -t30 * t33 + t37;
t1 = [MDP(1) + t28 * MDP(4) + (t38 * pkin(5) ^ 2 + t24 ^ 2) * MDP(14) + (t23 ^ 2 + t25 ^ 2 + t26 ^ 2) * MDP(18) + 0.2e1 * t38 * MDP(11) * pkin(5) + 0.2e1 * (t24 * MDP(12) + MDP(15) * t26 - t23 * MDP(17) + pkin(1) * MDP(9)) * t33 + 0.2e1 * (-pkin(1) * MDP(10) - t24 * MDP(13) + MDP(15) * t25 - t23 * MDP(16) + t33 * MDP(5)) * t32; t32 * MDP(6) + t33 * MDP(7) + (-t32 * pkin(2) + t27) * MDP(11) + (-t30 * t32 + t27) * MDP(15) + t26 * MDP(16) - t25 * MDP(17) + (t26 * qJ(3) - t25 * t30) * MDP(18) + ((qJ(3) * MDP(14) - MDP(10) + MDP(13)) * t33 + (-MDP(9) + t36) * t32) * pkin(5); MDP(8) - 0.2e1 * pkin(2) * MDP(12) + (pkin(2) ^ 2 + t34) * MDP(14) + 0.2e1 * t30 * MDP(17) + (t30 ^ 2 + t34) * MDP(18) + 0.2e1 * (MDP(13) + MDP(16)) * qJ(3); t25 * MDP(18) + (MDP(14) * pkin(5) + MDP(11) + MDP(15)) * t32; -t30 * MDP(18) - MDP(17) + t36; MDP(14) + MDP(18); t33 * MDP(15) + t26 * MDP(18); MDP(18) * qJ(3) + MDP(16); 0; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
