% Calculate Gravitation load on the joints for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4RRPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:35
% EndTime: 2019-07-18 18:16:35
% DurationCPUTime: 0.09s
% Computational Cost: add. (102->25), mult. (91->35), div. (0->0), fcn. (80->6), ass. (0->15)
t41 = qJ(1) + qJ(2);
t39 = sin(t41);
t40 = cos(t41);
t49 = sin(qJ(4));
t50 = cos(qJ(4));
t25 = -t39 * t49 - t40 * t50;
t26 = -t39 * t50 + t40 * t49;
t51 = (g(1) * t26 - g(2) * t25) * MDP(11) - (g(1) * t25 + g(2) * t26) * MDP(12);
t48 = t40 * pkin(2) + t39 * qJ(3);
t47 = -t39 * pkin(2) + t40 * qJ(3);
t31 = g(1) * t39 - g(2) * t40;
t44 = -t51 + (MDP(6) - MDP(8)) * (g(1) * t40 + g(2) * t39) + (MDP(5) + MDP(7)) * t31;
t43 = cos(qJ(1));
t42 = sin(qJ(1));
t1 = [(g(1) * t42 - g(2) * t43) * MDP(2) + (g(1) * t43 + g(2) * t42) * MDP(3) + (-g(1) * (-t42 * pkin(1) + t47) - g(2) * (t43 * pkin(1) + t48)) * MDP(9) + t44; (-g(1) * t47 - g(2) * t48) * MDP(9) + t44; -t31 * MDP(9); t51;];
taug  = t1;
