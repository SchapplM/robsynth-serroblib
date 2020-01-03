% Calculate Gravitation load on the joints for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPPRR12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:21
% EndTime: 2019-12-31 18:07:22
% DurationCPUTime: 0.15s
% Computational Cost: add. (85->40), mult. (130->54), div. (0->0), fcn. (112->8), ass. (0->23)
t41 = pkin(8) + qJ(4);
t36 = cos(t41);
t57 = t36 * MDP(17);
t44 = sin(qJ(5));
t46 = cos(qJ(5));
t56 = t46 * MDP(23) - t44 * MDP(24) + MDP(16);
t55 = g(3) * t36;
t45 = sin(qJ(1));
t54 = t45 * t44;
t53 = t45 * t46;
t47 = cos(qJ(1));
t52 = t47 * t44;
t51 = t47 * t46;
t50 = t47 * pkin(1) + t45 * qJ(2);
t34 = g(1) * t47 + g(2) * t45;
t33 = g(1) * t45 - g(2) * t47;
t38 = t47 * qJ(2);
t35 = sin(t41);
t32 = t35 * t51 - t54;
t31 = t35 * t52 + t53;
t30 = t35 * t53 + t52;
t29 = -t35 * t54 + t51;
t1 = [(-g(1) * (-pkin(1) * t45 + t38) - g(2) * t50) * MDP(6) + (-g(1) * (t38 + (-pkin(1) - qJ(3)) * t45) - g(2) * (qJ(3) * t47 + t50)) * MDP(10) + (-g(1) * t32 - g(2) * t30) * MDP(23) + (g(1) * t31 - g(2) * t29) * MDP(24) + (MDP(2) - MDP(4) + MDP(9)) * t33 + (-t35 * MDP(16) - t57 - MDP(7) * sin(pkin(8)) - MDP(8) * cos(pkin(8)) + MDP(3) - MDP(5)) * t34; (-MDP(10) - MDP(6)) * t33; -t34 * MDP(10); (t35 * t56 + t57) * g(3) + (MDP(17) * t35 - t36 * t56) * t33; (-g(1) * t29 - g(2) * t31 + t44 * t55) * MDP(23) + (g(1) * t30 - g(2) * t32 + t46 * t55) * MDP(24);];
taug = t1;
