% Calculate Gravitation load on the joints for
% S4RRPP2
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
%   see S4RRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4RRPP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:34:05
% EndTime: 2019-03-08 18:34:06
% DurationCPUTime: 0.05s
% Computational Cost: add. (101->30), mult. (85->34), div. (0->0), fcn. (56->4), ass. (0->15)
t36 = sin(qJ(1));
t43 = t36 * pkin(1);
t35 = qJ(1) + qJ(2);
t32 = sin(t35);
t33 = cos(t35);
t42 = t33 * pkin(2) + t32 * qJ(3);
t37 = cos(qJ(1));
t41 = t37 * pkin(1) + t42;
t28 = t33 * qJ(3);
t40 = -t32 * pkin(2) + t28;
t39 = t28 + (-pkin(2) - pkin(3)) * t32;
t25 = g(1) * t32 - g(2) * t33;
t38 = (-MDP(11) + MDP(6) - MDP(8)) * (g(1) * t33 + g(2) * t32) + (MDP(10) + MDP(5) + MDP(7)) * t25;
t29 = t33 * pkin(3);
t1 = [(g(1) * t36 - g(2) * t37) * MDP(2) + (g(1) * t37 + g(2) * t36) * MDP(3) + (-g(1) * (t40 - t43) - g(2) * t41) * MDP(9) + (-g(1) * (t39 - t43) - g(2) * (t29 + t41)) * MDP(12) + t38; (-g(1) * t40 - g(2) * t42) * MDP(9) + (-g(1) * t39 - g(2) * (t29 + t42)) * MDP(12) + t38; (-MDP(12) - MDP(9)) * t25; g(3) * MDP(12);];
taug  = t1;
