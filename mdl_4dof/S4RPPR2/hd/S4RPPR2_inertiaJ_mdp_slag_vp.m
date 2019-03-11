% Calculate joint inertia matrix for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4RPPR2_inertiaJ_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:28:32
% EndTime: 2019-03-08 18:28:33
% DurationCPUTime: 0.04s
% Computational Cost: add. (63->30), mult. (81->41), div. (0->0), fcn. (58->4), ass. (0->12)
t31 = -pkin(1) - pkin(2);
t24 = sin(pkin(6));
t25 = cos(pkin(6));
t20 = t24 * qJ(2) - t25 * t31;
t19 = -pkin(3) - t20;
t22 = t25 * qJ(2) + t24 * t31;
t26 = sin(qJ(4));
t27 = cos(qJ(4));
t30 = (-t27 * t19 + t26 * t22) * MDP(11);
t29 = (t26 * t19 + t27 * t22) * MDP(12);
t28 = -(t26 * t24 - t27 * t25) * MDP(11) - (t27 * t24 + t26 * t25) * MDP(12);
t1 = [MDP(1) + (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + (t20 ^ 2 + t22 ^ 2) * MDP(9) + MDP(10) + 0.2e1 * t30 + 0.2e1 * t29 + 0.2e1 * t20 * MDP(7) + 0.2e1 * t22 * MDP(8); -MDP(4) - pkin(1) * MDP(6) - t25 * MDP(7) + t24 * MDP(8) + (-t20 * t25 + t22 * t24) * MDP(9) - t28; MDP(6) + (t24 ^ 2 + t25 ^ 2) * MDP(9); 0; 0; MDP(9); -MDP(10) - t29 - t30; t28; 0; MDP(10);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
