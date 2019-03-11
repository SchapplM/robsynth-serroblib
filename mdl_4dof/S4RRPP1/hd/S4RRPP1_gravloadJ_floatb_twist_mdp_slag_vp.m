% Calculate Gravitation load on the joints for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RRPP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:33:10
% EndTime: 2019-03-08 18:33:10
% DurationCPUTime: 0.05s
% Computational Cost: add. (92->30), mult. (67->39), div. (0->0), fcn. (42->6), ass. (0->18)
t37 = qJ(1) + qJ(2);
t34 = sin(t37);
t45 = pkin(2) * t34;
t38 = sin(qJ(1));
t44 = t38 * pkin(1);
t33 = pkin(6) + t37;
t30 = sin(t33);
t31 = cos(t33);
t35 = cos(t37);
t32 = pkin(2) * t35;
t43 = t31 * pkin(3) + t30 * qJ(4) + t32;
t24 = g(1) * t30 - g(2) * t31;
t41 = g(1) * t34 - g(2) * t35;
t42 = t24 * MDP(8) + (-g(1) * t31 - g(2) * t30) * MDP(9) + t41 * MDP(5) + (g(1) * t35 + g(2) * t34) * MDP(6);
t40 = -t30 * pkin(3) + t31 * qJ(4) - t45;
t39 = cos(qJ(1));
t36 = t39 * pkin(1);
t1 = [(g(1) * t38 - g(2) * t39) * MDP(2) + (g(1) * t39 + g(2) * t38) * MDP(3) + (-g(1) * (-t44 - t45) - g(2) * (t32 + t36)) * MDP(7) + (-g(1) * (t40 - t44) - g(2) * (t36 + t43)) * MDP(10) + t42; t41 * pkin(2) * MDP(7) + (-g(1) * t40 - g(2) * t43) * MDP(10) + t42; (-MDP(10) - MDP(7)) * g(3); -t24 * MDP(10);];
taug  = t1;
