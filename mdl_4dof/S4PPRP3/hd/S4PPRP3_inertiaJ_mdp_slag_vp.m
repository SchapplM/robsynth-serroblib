% Calculate joint inertia matrix for
% S4PPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
% MDP [6x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP3_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP3_inertiaJ_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'S4PPRP3_inertiaJ_mdp_slag_vp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:14:40
% EndTime: 2019-03-08 18:14:40
% DurationCPUTime: 0.02s
% Computational Cost: add. (10->7), mult. (16->10), div. (0->0), fcn. (10->2), ass. (0->5)
t5 = sin(qJ(3));
t6 = cos(qJ(3));
t8 = MDP(2) + (t5 ^ 2 + t6 ^ 2) * MDP(6);
t7 = MDP(6) * pkin(3) + MDP(4);
t1 = [MDP(1) + t8; 0; t8; -t6 * MDP(5) - t7 * t5; -t5 * MDP(5) + t7 * t6; pkin(3) ^ 2 * MDP(6) + MDP(3); 0; 0; 0; MDP(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
