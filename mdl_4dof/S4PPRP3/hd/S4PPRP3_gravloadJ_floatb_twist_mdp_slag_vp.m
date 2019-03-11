% Calculate Gravitation load on the joints for
% S4PPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
% MDP [6x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PPRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'S4PPRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:14:40
% EndTime: 2019-03-08 18:14:40
% DurationCPUTime: 0.02s
% Computational Cost: add. (8->6), mult. (17->10), div. (0->0), fcn. (6->2), ass. (0->4)
t7 = -MDP(2) - MDP(6);
t5 = cos(qJ(3));
t4 = sin(qJ(3));
t1 = [(-MDP(1) + t7) * g(2); t7 * g(1); (g(1) * t4 + g(2) * t5) * MDP(5) + (MDP(6) * pkin(3) + MDP(4)) * (-g(1) * t5 + g(2) * t4); g(3) * MDP(6);];
taug  = t1;
