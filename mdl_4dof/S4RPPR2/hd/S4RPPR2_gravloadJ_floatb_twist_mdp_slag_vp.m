% Calculate Gravitation load on the joints for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4RPPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:28:32
% EndTime: 2019-03-08 18:28:32
% DurationCPUTime: 0.07s
% Computational Cost: add. (59->28), mult. (77->40), div. (0->0), fcn. (70->6), ass. (0->16)
t39 = sin(qJ(1));
t40 = cos(qJ(1));
t45 = pkin(6) + qJ(4);
t43 = sin(t45);
t44 = cos(t45);
t24 = -t39 * t43 - t40 * t44;
t25 = -t39 * t44 + t40 * t43;
t47 = (g(1) * t25 - g(2) * t24) * MDP(11) - (g(1) * t24 + g(2) * t25) * MDP(12);
t46 = t40 * pkin(1) + t39 * qJ(2);
t38 = cos(pkin(6));
t37 = sin(pkin(6));
t34 = t40 * qJ(2);
t28 = g(1) * t39 - g(2) * t40;
t27 = t39 * t37 + t40 * t38;
t26 = t40 * t37 - t39 * t38;
t1 = [(-g(1) * (-t39 * pkin(1) + t34) - g(2) * t46) * MDP(6) + (-g(1) * t26 - g(2) * t27) * MDP(7) + (-g(1) * t27 + g(2) * t26) * MDP(8) + (-g(1) * (t34 + (-pkin(1) - pkin(2)) * t39) - g(2) * (t40 * pkin(2) + t46)) * MDP(9) + (MDP(3) - MDP(5)) * (g(1) * t40 + g(2) * t39) + (MDP(2) + MDP(4)) * t28 - t47; (-MDP(6) - MDP(9)) * t28; g(3) * MDP(9); t47;];
taug  = t1;
