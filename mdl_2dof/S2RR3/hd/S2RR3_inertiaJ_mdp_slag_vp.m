% Calculate joint inertia matrix for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% MDP [6x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S2RR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S2RR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_inertiaJ_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_inertiaJ_mdp_slag_vp: pkin has to be [3x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'S2RR3_inertiaJ_mdp_slag_vp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:27
% EndTime: 2020-06-19 09:14:28
% DurationCPUTime: 0.03s
% Computational Cost: add. (5->4), mult. (10->4), div. (0->0), fcn. (4->2), ass. (0->2)
t7 = (MDP(5) * cos(qJ(2)) - MDP(6) * sin(qJ(2))) * pkin(1);
t1 = [MDP(1) + MDP(4) + 0.2e1 * t7; MDP(4) + t7; MDP(4);];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t1(1), t1(2); t1(2), t1(3);];
Mq = res;
