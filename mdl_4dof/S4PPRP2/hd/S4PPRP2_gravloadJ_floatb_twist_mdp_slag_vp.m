% Calculate Gravitation load on the joints for
% S4PPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta2]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PPRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PPRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:13:21
% EndTime: 2019-03-08 18:13:21
% DurationCPUTime: 0.03s
% Computational Cost: add. (29->13), mult. (27->16), div. (0->0), fcn. (14->2), ass. (0->6)
t11 = -MDP(2) - MDP(8);
t10 = pkin(5) + qJ(3);
t9 = cos(t10);
t8 = sin(t10);
t6 = g(1) * t8 - g(2) * t9;
t1 = [(-MDP(1) + t11) * g(2); t11 * g(3); (-g(1) * (-t8 * pkin(3) + t9 * qJ(4)) - g(2) * (t9 * pkin(3) + t8 * qJ(4))) * MDP(8) + (MDP(5) - MDP(7)) * (g(1) * t9 + g(2) * t8) + (MDP(4) + MDP(6)) * t6; -t6 * MDP(8);];
taug  = t1;
