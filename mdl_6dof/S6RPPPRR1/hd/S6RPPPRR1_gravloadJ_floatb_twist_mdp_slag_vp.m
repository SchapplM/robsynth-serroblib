% Calculate Gravitation load on the joints for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPPRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPPPRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:30:14
% EndTime: 2019-03-09 01:30:14
% DurationCPUTime: 0.16s
% Computational Cost: add. (117->45), mult. (135->64), div. (0->0), fcn. (112->8), ass. (0->22)
t49 = cos(qJ(5));
t59 = g(3) * t49;
t45 = sin(qJ(6));
t46 = sin(qJ(5));
t58 = t45 * t46;
t48 = cos(qJ(6));
t57 = t46 * t48;
t56 = -MDP(10) - MDP(7);
t44 = qJ(1) + pkin(9);
t41 = sin(t44);
t42 = cos(t44);
t50 = cos(qJ(1));
t55 = t50 * pkin(1) + t42 * pkin(2) + t41 * qJ(3);
t47 = sin(qJ(1));
t54 = -t47 * pkin(1) + t42 * qJ(3);
t53 = g(1) * t42 + g(2) * t41;
t35 = g(1) * t41 - g(2) * t42;
t34 = -t41 * t45 + t42 * t57;
t33 = -t41 * t48 - t42 * t58;
t32 = -t41 * t57 - t42 * t45;
t31 = t41 * t58 - t42 * t48;
t1 = [(g(1) * t50 + g(2) * t47) * MDP(3) + (-g(1) * (-t41 * pkin(2) + t54) - g(2) * t55) * MDP(7) + (-g(1) * ((-pkin(2) - qJ(4)) * t41 + t54) - g(2) * (t42 * qJ(4) + t55)) * MDP(10) + (-g(1) * t32 - g(2) * t34) * MDP(23) + (-g(1) * t31 - g(2) * t33) * MDP(24) - (MDP(6) + MDP(8)) * t53 + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t47 - g(2) * t50) + (MDP(16) * t46 + MDP(17) * t49 - MDP(5) + MDP(9)) * t35; (-MDP(4) + t56) * g(3); t56 * t35; -t53 * MDP(10); (t53 * t46 + t59) * MDP(17) + (-MDP(23) * t48 + MDP(24) * t45 - MDP(16)) * (-g(3) * t46 + t53 * t49); (-g(1) * t33 + g(2) * t31 + t45 * t59) * MDP(23) + (g(1) * t34 - g(2) * t32 + t48 * t59) * MDP(24);];
taug  = t1;
