% Calculate joint inertia matrix for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RPRR5_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:40
% EndTime: 2019-12-31 16:51:40
% DurationCPUTime: 0.06s
% Computational Cost: add. (64->31), mult. (97->39), div. (0->0), fcn. (61->4), ass. (0->21)
t26 = sin(qJ(3));
t28 = cos(qJ(3));
t25 = sin(qJ(4));
t27 = cos(qJ(4));
t36 = t27 * MDP(15);
t33 = -t25 * MDP(16) + t36;
t43 = (MDP(8) + t33) * t28 - t26 * MDP(9);
t29 = -pkin(1) - pkin(2);
t20 = t26 * qJ(2) - t28 * t29;
t18 = pkin(3) + t20;
t42 = pkin(3) + t18;
t41 = t20 * MDP(8);
t21 = t28 * qJ(2) + t26 * t29;
t40 = t21 * MDP(9);
t38 = t25 ^ 2 * MDP(10) + MDP(7);
t37 = t27 * MDP(11);
t35 = 0.2e1 * t25 * t37 + t38;
t34 = -t25 * MDP(12) - t27 * MDP(13);
t32 = -MDP(15) * t25 - MDP(16) * t27;
t30 = 0.2e1 * t33;
t1 = [MDP(1) + (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + t18 * t30 + 0.2e1 * t41 + 0.2e1 * t40 + t35; -pkin(1) * MDP(6) - MDP(4) - t43; MDP(6); -t41 - t40 - t42 * t36 + (t42 * MDP(16) - 0.2e1 * t37) * t25 - t38; t43; pkin(3) * t30 + t35; t32 * (-pkin(6) + t21) + t34; t32 * t26; t32 * pkin(6) - t34; MDP(14);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
