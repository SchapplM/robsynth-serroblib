% Calculate joint inertia matrix for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RPPR1_inertiaJ_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:27:31
% EndTime: 2019-03-08 18:27:32
% DurationCPUTime: 0.03s
% Computational Cost: add. (37->20), mult. (50->25), div. (0->0), fcn. (27->4), ass. (0->11)
t19 = cos(pkin(6));
t16 = t19 * pkin(1) + pkin(2);
t14 = -pkin(3) - t16;
t18 = sin(pkin(6));
t15 = t18 * pkin(1) + qJ(3);
t20 = sin(qJ(4));
t21 = cos(qJ(4));
t25 = (-t21 * t14 + t20 * t15) * MDP(9);
t24 = (t20 * t14 + t21 * t15) * MDP(10);
t23 = -t20 * MDP(10) + t21 * MDP(9);
t1 = [MDP(1) + (t15 ^ 2 + t16 ^ 2) * MDP(7) + MDP(8) + (t18 ^ 2 + t19 ^ 2) * MDP(4) * pkin(1) ^ 2 + 0.2e1 * t24 + 0.2e1 * t16 * MDP(5) + 0.2e1 * t15 * MDP(6) + 0.2e1 * t25; 0; MDP(4) + MDP(7); -t16 * MDP(7) - MDP(5) - t23; 0; MDP(7); -MDP(8) - t24 - t25; 0; t23; MDP(8);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
