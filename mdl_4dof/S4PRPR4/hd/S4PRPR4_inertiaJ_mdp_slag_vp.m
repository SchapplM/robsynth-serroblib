% Calculate joint inertia matrix for
% S4PRPR4
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
%   see S4PRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4PRPR4_inertiaJ_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:00
% EndTime: 2019-12-31 16:22:00
% DurationCPUTime: 0.02s
% Computational Cost: add. (17->16), mult. (29->20), div. (0->0), fcn. (13->2), ass. (0->7)
t10 = (pkin(2) * MDP(7));
t5 = sin(qJ(4));
t9 = t5 * MDP(13);
t6 = cos(qJ(4));
t8 = t6 * MDP(14);
t7 = -pkin(2) - pkin(5);
t1 = [MDP(1) + MDP(7); 0; MDP(2) + (MDP(8) * t6 - 0.2e1 * t5 * MDP(9)) * t6 + ((-2 * MDP(5) + t10) * pkin(2)) + (MDP(7) * qJ(3) + (2 * MDP(6)) + 0.2e1 * t8 + 0.2e1 * t9) * qJ(3); 0; MDP(5) - t10; MDP(7); -t8 - t9; (MDP(13) * t7 + MDP(10)) * t6 + (-MDP(14) * t7 - MDP(11)) * t5; t6 * MDP(13) - t5 * MDP(14); MDP(12);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
