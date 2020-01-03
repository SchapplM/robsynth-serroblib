% Calculate joint inertia matrix for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4PRRR3_inertiaJ_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:39
% EndTime: 2019-12-31 16:31:39
% DurationCPUTime: 0.04s
% Computational Cost: add. (32->17), mult. (63->23), div. (0->0), fcn. (37->4), ass. (0->11)
t19 = sin(qJ(4));
t21 = cos(qJ(4));
t25 = t21 * MDP(13) - t19 * MDP(14);
t29 = t19 * MDP(10) + t21 * MDP(11);
t26 = MDP(5) + (MDP(8) * t19 + 0.2e1 * MDP(9) * t21) * t19;
t24 = -MDP(13) * t19 - MDP(14) * t21;
t20 = sin(qJ(3));
t22 = cos(qJ(3));
t23 = (t22 * MDP(6) - t20 * MDP(7)) * pkin(2);
t15 = -t22 * pkin(2) - pkin(3);
t1 = [MDP(1); 0; -0.2e1 * t15 * t25 + MDP(2) + 0.2e1 * t23 + t26; 0; t23 + t26 + t25 * (pkin(3) - t15); 0.2e1 * pkin(3) * t25 + t26; t25; t24 * (t20 * pkin(2) + pkin(6)) + t29; t24 * pkin(6) + t29; MDP(12);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
