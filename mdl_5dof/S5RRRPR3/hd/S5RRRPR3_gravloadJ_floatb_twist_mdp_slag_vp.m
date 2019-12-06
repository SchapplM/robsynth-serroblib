% Calculate Gravitation load on the joints for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:42:58
% EndTime: 2019-12-05 18:42:58
% DurationCPUTime: 0.07s
% Computational Cost: add. (160->33), mult. (132->43), div. (0->0), fcn. (99->8), ass. (0->19)
t44 = qJ(1) + qJ(2);
t42 = sin(t44);
t43 = cos(t44);
t35 = g(2) * t42 - g(3) * t43;
t41 = qJ(3) + pkin(9) + qJ(5);
t38 = sin(t41);
t39 = cos(t41);
t54 = (-g(1) * t39 - t35 * t38) * MDP(21) + (g(1) * t38 - t35 * t39) * MDP(22);
t48 = cos(qJ(3));
t40 = t48 * pkin(3) + pkin(2);
t45 = -qJ(4) - pkin(7);
t53 = -t43 * t40 + t42 * t45;
t36 = g(2) * t43 + g(3) * t42;
t52 = -t42 * t40 - t43 * t45;
t46 = sin(qJ(3));
t51 = (MDP(14) - MDP(6)) * t35 + (t48 * MDP(12) - t46 * MDP(13) + t39 * MDP(21) - t38 * MDP(22) + MDP(5)) * t36;
t49 = cos(qJ(1));
t47 = sin(qJ(1));
t1 = [(g(2) * t49 + g(3) * t47) * MDP(2) + (-g(2) * t47 + g(3) * t49) * MDP(3) + (-g(2) * (-t49 * pkin(1) + t53) - g(3) * (-t47 * pkin(1) + t52)) * MDP(15) + t51; (-g(2) * t53 - g(3) * t52) * MDP(15) + t51; (g(1) * t46 - t35 * t48) * MDP(13) + t54 + (MDP(15) * pkin(3) + MDP(12)) * (-g(1) * t48 - t35 * t46); -t36 * MDP(15); t54;];
taug = t1;
