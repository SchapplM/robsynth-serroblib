% Calculate joint inertia matrix for
% S4PPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PPRR2_inertiaJ_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:00
% EndTime: 2019-03-08 18:17:00
% DurationCPUTime: 0.03s
% Computational Cost: add. (25->13), mult. (47->19), div. (0->0), fcn. (54->6), ass. (0->11)
t17 = sin(pkin(6));
t18 = cos(pkin(6));
t20 = sin(qJ(3));
t22 = cos(qJ(3));
t15 = -t20 * t17 + t22 * t18;
t16 = t22 * t17 + t20 * t18;
t19 = sin(qJ(4));
t21 = cos(qJ(4));
t24 = (t21 * t15 - t19 * t16) * MDP(7) + (-t19 * t15 - t21 * t16) * MDP(8);
t23 = (MDP(7) * t21 - MDP(8) * t19) * pkin(3);
t1 = [MDP(1) + (t17 ^ 2 + t18 ^ 2) * MDP(2); 0; MDP(2); t15 * MDP(4) - t16 * MDP(5) + t24; 0; MDP(3) + MDP(6) + 0.2e1 * t23; t24; 0; MDP(6) + t23; MDP(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
