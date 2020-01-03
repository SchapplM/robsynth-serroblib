% Calculate joint inertia matrix for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4PRPR7_inertiaJ_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:55
% EndTime: 2019-12-31 16:25:55
% DurationCPUTime: 0.03s
% Computational Cost: add. (26->23), mult. (46->26), div. (0->0), fcn. (28->4), ass. (0->11)
t17 = (MDP(7) * pkin(2));
t7 = sin(qJ(4));
t16 = t7 * MDP(13);
t9 = cos(qJ(4));
t15 = t9 * MDP(14);
t14 = qJ(3) * MDP(7);
t13 = MDP(5) - t17;
t12 = t9 * MDP(13) - t7 * MDP(14);
t10 = cos(qJ(2));
t8 = sin(qJ(2));
t1 = [MDP(1) + (t10 ^ 2 + t8 ^ 2) * MDP(7); (MDP(3) - t13) * t10 + (-MDP(4) + MDP(6) + t14 + t15 + t16) * t8; MDP(2) + (MDP(8) * t9 - 0.2e1 * t7 * MDP(9)) * t9 + ((-2 * MDP(5) + t17) * pkin(2)) + (0.2e1 * MDP(6) + t14 + 0.2e1 * t15 + 0.2e1 * t16) * qJ(3); -t10 * MDP(7); t13; MDP(7); -t12 * t10; t9 * MDP(10) - t7 * MDP(11) + t12 * (-pkin(2) - pkin(5)); t12; MDP(12);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
