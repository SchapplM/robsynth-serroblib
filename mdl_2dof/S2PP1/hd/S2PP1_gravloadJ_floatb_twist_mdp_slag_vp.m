% Calculate Gravitation load on the joints for
% S2PP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2]';
% MDP [2x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S2PP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-03-03 18:41
% Revision: 33b345ae0dd6ec4aa15499ab3d43edbbded0bea5 (2021-02-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S2PP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(1,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2PP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2PP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [1x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [2 1]), ...
  'S2PP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [2x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-03-03 18:41:23
% EndTime: 2021-03-03 18:41:24
% DurationCPUTime: 0.02s
% Computational Cost: add. (2->2), mult. (3->2), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [(-MDP(1) - MDP(2)) * g(1); -g(2) * MDP(2);];
taug = t1;
