% Calculate Gravitation load on the joints for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:51
% EndTime: 2019-12-05 16:15:51
% DurationCPUTime: 0.06s
% Computational Cost: add. (151->28), mult. (97->36), div. (0->0), fcn. (70->8), ass. (0->15)
t43 = pkin(8) + qJ(2);
t41 = qJ(3) + t43;
t35 = sin(t41);
t36 = cos(t41);
t48 = t36 * pkin(3) + t35 * qJ(4);
t47 = -t35 * pkin(3) + t36 * qJ(4);
t31 = g(1) * t36 + g(2) * t35;
t30 = g(1) * t35 - g(2) * t36;
t42 = pkin(9) + qJ(5);
t37 = sin(t42);
t39 = cos(t42);
t46 = (-MDP(10) + MDP(7)) * t31 + (t39 * MDP(17) - t37 * MDP(18) + MDP(8) * cos(pkin(9)) - MDP(9) * sin(pkin(9)) + MDP(6)) * t30;
t40 = cos(t43);
t38 = sin(t43);
t1 = [(-MDP(1) - MDP(11)) * g(3); (g(1) * t38 - g(2) * t40) * MDP(3) + (g(1) * t40 + g(2) * t38) * MDP(4) + (-g(1) * (-pkin(2) * t38 + t47) - g(2) * (pkin(2) * t40 + t48)) * MDP(11) + t46; (-g(1) * t47 - g(2) * t48) * MDP(11) + t46; -t30 * MDP(11); (-g(3) * t39 + t31 * t37) * MDP(17) + (g(3) * t37 + t31 * t39) * MDP(18);];
taug = t1;
