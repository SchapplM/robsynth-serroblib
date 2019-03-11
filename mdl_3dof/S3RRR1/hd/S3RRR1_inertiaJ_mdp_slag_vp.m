% Calculate joint inertia matrix for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% MDP [9x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3RRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [3x3]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S3RRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_inertiaJ_mdp_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [9 1]), ...
  'S3RRR1_inertiaJ_mdp_slag_vp: MDP has to be [9x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:08:08
% EndTime: 2019-03-08 18:08:08
% DurationCPUTime: 0.04s
% Computational Cost: add. (32->22), mult. (56->25), div. (0->0), fcn. (34->4), ass. (0->14)
t19 = sin(qJ(2));
t28 = pkin(1) * t19;
t21 = cos(qJ(2));
t17 = t21 * pkin(1) + pkin(2);
t20 = cos(qJ(3));
t16 = t20 * t17;
t18 = sin(qJ(3));
t27 = (-t18 * t28 + t16) * MDP(8);
t26 = (-t18 * t17 - t20 * t28) * MDP(9);
t25 = t18 * MDP(9);
t24 = t21 * MDP(5);
t23 = MDP(4) + MDP(7);
t22 = (t20 * MDP(8) - t25) * pkin(2);
t1 = [MDP(1) + 0.2e1 * (-t19 * MDP(6) + t24) * pkin(1) + 0.2e1 * t27 + 0.2e1 * t26 + t23; (t20 * pkin(2) + t16) * MDP(8) + (-pkin(2) - t17) * t25 + (t24 + (-MDP(8) * t18 - MDP(9) * t20 - MDP(6)) * t19) * pkin(1) + t23; 0.2e1 * t22 + t23; MDP(7) + t26 + t27; MDP(7) + t22; MDP(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t1(1) t1(2) t1(4); t1(2) t1(3) t1(5); t1(4) t1(5) t1(6);];
Mq  = res;
