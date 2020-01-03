% Calculate Gravitation load on the joints for
% S5RPRRR12
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
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RPRRR12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:13:09
% EndTime: 2019-12-31 19:13:09
% DurationCPUTime: 0.13s
% Computational Cost: add. (103->37), mult. (154->54), div. (0->0), fcn. (136->8), ass. (0->21)
t41 = qJ(3) + qJ(4);
t39 = cos(t41);
t55 = g(3) * t39;
t42 = sin(qJ(5));
t44 = sin(qJ(1));
t54 = t44 * t42;
t45 = cos(qJ(5));
t53 = t44 * t45;
t47 = cos(qJ(1));
t52 = t47 * t42;
t51 = t47 * t45;
t36 = g(1) * t44 - g(2) * t47;
t38 = sin(t41);
t50 = (t36 * t38 + t55) * MDP(20) + (-t45 * MDP(26) + t42 * MDP(27) - MDP(19)) * (-g(3) * t38 + t36 * t39);
t46 = cos(qJ(3));
t43 = sin(qJ(3));
t35 = t38 * t51 - t54;
t34 = t38 * t52 + t53;
t33 = t38 * t53 + t52;
t32 = -t38 * t54 + t51;
t1 = [(-g(1) * (-t44 * pkin(1) + t47 * qJ(2)) - g(2) * (t47 * pkin(1) + t44 * qJ(2))) * MDP(6) + (-g(1) * t35 - g(2) * t33) * MDP(26) + (g(1) * t34 - g(2) * t32) * MDP(27) + (MDP(2) - MDP(4)) * t36 + (-t43 * MDP(12) - t46 * MDP(13) - MDP(19) * t38 - MDP(20) * t39 + MDP(3) - MDP(5)) * (g(1) * t47 + g(2) * t44); -t36 * MDP(6); (g(3) * t43 - t36 * t46) * MDP(12) + (g(3) * t46 + t36 * t43) * MDP(13) + t50; t50; (-g(1) * t32 - g(2) * t34 + t42 * t55) * MDP(26) + (g(1) * t33 - g(2) * t35 + t45 * t55) * MDP(27);];
taug = t1;
