% Calculate Gravitation load on the joints for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:54:42
% EndTime: 2019-12-05 15:54:43
% DurationCPUTime: 0.18s
% Computational Cost: add. (131->38), mult. (160->57), div. (0->0), fcn. (150->10), ass. (0->17)
t37 = sin(pkin(8));
t39 = cos(pkin(8));
t44 = g(1) * t39 + g(2) * t37;
t40 = sin(qJ(2));
t41 = cos(qJ(2));
t28 = -g(3) * t41 + t44 * t40;
t49 = g(3) * t40;
t47 = t37 * t41;
t46 = t39 * t41;
t35 = pkin(9) + qJ(4);
t34 = qJ(5) + t35;
t30 = sin(t34);
t31 = cos(t34);
t45 = (-g(1) * (-t30 * t46 + t37 * t31) - g(2) * (-t30 * t47 - t39 * t31) + t30 * t49) * MDP(21) + (-g(1) * (-t37 * t30 - t31 * t46) - g(2) * (t39 * t30 - t31 * t47) + t31 * t49) * MDP(22);
t33 = cos(t35);
t32 = sin(t35);
t1 = [(-MDP(1) - MDP(8)) * g(3); (-g(3) * (t41 * pkin(2) + t40 * qJ(3)) + t44 * (pkin(2) * t40 - qJ(3) * t41)) * MDP(8) + (MDP(4) - MDP(7)) * (t44 * t41 + t49) + (t33 * MDP(14) - t32 * MDP(15) + t31 * MDP(21) - t30 * MDP(22) + cos(pkin(9)) * MDP(5) - sin(pkin(9)) * MDP(6) + MDP(3)) * t28; -t28 * MDP(8); (-g(1) * (-t32 * t46 + t37 * t33) - g(2) * (-t32 * t47 - t39 * t33) + t32 * t49) * MDP(14) + (-g(1) * (-t37 * t32 - t33 * t46) - g(2) * (t39 * t32 - t33 * t47) + t33 * t49) * MDP(15) + t45; t45;];
taug = t1;
