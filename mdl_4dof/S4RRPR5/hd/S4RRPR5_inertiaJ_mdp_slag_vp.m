% Calculate joint inertia matrix for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RRPR5_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:33
% EndTime: 2019-12-31 17:03:33
% DurationCPUTime: 0.06s
% Computational Cost: add. (59->38), mult. (97->42), div. (0->0), fcn. (46->4), ass. (0->19)
t40 = -2 * MDP(7);
t39 = 2 * MDP(8);
t38 = (MDP(9) * pkin(2));
t28 = cos(qJ(2));
t22 = -t28 * pkin(1) - pkin(2);
t37 = t22 * MDP(9);
t25 = sin(qJ(4));
t35 = t25 * MDP(15);
t27 = cos(qJ(4));
t34 = t27 * MDP(15);
t33 = t27 * MDP(16);
t32 = MDP(4) + (MDP(10) * t27 - 0.2e1 * MDP(11) * t25) * t27;
t31 = -t25 * MDP(16) + t34;
t30 = t39 + 0.2e1 * t33 + 0.2e1 * t35;
t29 = -pkin(2) - pkin(6);
t26 = sin(qJ(2));
t24 = t27 * MDP(12);
t20 = t26 * pkin(1) + qJ(3);
t1 = [MDP(1) + ((2 * MDP(7)) + t37) * t22 + 0.2e1 * (t28 * MDP(5) - t26 * MDP(6)) * pkin(1) + (t20 * MDP(9) + t30) * t20 + t32; pkin(2) * t40 + qJ(3) * t39 + (-t22 * pkin(2) + t20 * qJ(3)) * MDP(9) + ((MDP(5) - MDP(7)) * t28 + (-MDP(6) + MDP(8)) * t26) * pkin(1) + t32 + (t33 + t35) * (qJ(3) + t20); (t40 + t38) * pkin(2) + (MDP(9) * qJ(3) + t30) * qJ(3) + t32; MDP(7) + t37; MDP(7) - t38; MDP(9); -t25 * MDP(13) + t24 + t31 * (-pkin(6) + t22); t29 * t34 + t24 + (-MDP(16) * t29 - MDP(13)) * t25; t31; MDP(14);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
