% Calculate Gravitation load on the joints for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRRR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:06
% EndTime: 2019-12-31 19:06:06
% DurationCPUTime: 0.11s
% Computational Cost: add. (107->28), mult. (182->40), div. (0->0), fcn. (192->8), ass. (0->16)
t46 = sin(qJ(3));
t47 = sin(qJ(1));
t48 = cos(qJ(3));
t49 = cos(qJ(1));
t27 = -t47 * t46 - t49 * t48;
t28 = t49 * t46 - t47 * t48;
t36 = qJ(4) + qJ(5);
t34 = sin(t36);
t35 = cos(t36);
t37 = sin(qJ(4));
t38 = cos(qJ(4));
t39 = g(1) * t27 + g(2) * t28;
t52 = (t38 * MDP(15) - t37 * MDP(16) + MDP(22) * t35 - MDP(23) * t34 + MDP(8)) * (g(1) * t28 - g(2) * t27) - t39 * MDP(9);
t45 = (g(3) * t35 - t39 * t34) * MDP(22) + (-g(3) * t34 - t39 * t35) * MDP(23);
t29 = g(1) * t47 - g(2) * t49;
t1 = [(-g(1) * (-t47 * pkin(1) + t49 * qJ(2)) - g(2) * (t49 * pkin(1) + t47 * qJ(2))) * MDP(6) + (MDP(3) - MDP(5)) * (g(1) * t49 + g(2) * t47) + (MDP(2) + MDP(4)) * t29 - t52; -t29 * MDP(6); t52; (g(3) * t38 - t39 * t37) * MDP(15) + (-g(3) * t37 - t39 * t38) * MDP(16) + t45; t45;];
taug = t1;
