% Calculate joint inertia matrix for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PPRR1_inertiaJ_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:15:49
% EndTime: 2019-03-08 18:15:49
% DurationCPUTime: 0.02s
% Computational Cost: add. (14->10), mult. (24->12), div. (0->0), fcn. (22->4), ass. (0->7)
t11 = sin(qJ(4));
t12 = sin(qJ(3));
t13 = cos(qJ(4));
t14 = cos(qJ(3));
t16 = (-t11 * t12 + t13 * t14) * MDP(7) + (-t11 * t14 - t13 * t12) * MDP(8);
t15 = (MDP(7) * t13 - MDP(8) * t11) * pkin(3);
t1 = [MDP(1) + MDP(2); 0; MDP(2); 0; t14 * MDP(4) - t12 * MDP(5) + t16; MDP(3) + MDP(6) + 0.2e1 * t15; 0; t16; MDP(6) + t15; MDP(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
