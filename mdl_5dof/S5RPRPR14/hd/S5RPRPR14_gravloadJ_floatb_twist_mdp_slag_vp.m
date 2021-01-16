% Calculate Gravitation load on the joints for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR14_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:17
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR14_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR14_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:17:04
% EndTime: 2021-01-15 12:17:06
% DurationCPUTime: 0.17s
% Computational Cost: add. (98->45), mult. (151->62), div. (0->0), fcn. (127->8), ass. (0->25)
t40 = sin(qJ(5));
t43 = cos(qJ(5));
t55 = t43 * MDP(23) - t40 * MDP(24) + MDP(14);
t39 = qJ(3) + pkin(8);
t36 = cos(t39);
t44 = cos(qJ(3));
t54 = t44 * MDP(13) + MDP(15) * t36;
t53 = g(3) * t36;
t42 = sin(qJ(1));
t52 = t42 * t40;
t51 = t42 * t43;
t45 = cos(qJ(1));
t50 = t45 * t40;
t49 = t45 * t43;
t33 = g(1) * t45 + g(2) * t42;
t32 = g(1) * t42 - g(2) * t45;
t41 = sin(qJ(3));
t38 = pkin(1) + pkin(6) + qJ(4);
t35 = sin(t39);
t34 = t41 * pkin(3) + qJ(2);
t31 = t35 * t49 - t52;
t30 = t35 * t50 + t51;
t29 = t35 * t51 + t50;
t28 = -t35 * t52 + t49;
t1 = [(-g(1) * (-t42 * pkin(1) + t45 * qJ(2)) - g(2) * (t45 * pkin(1) + t42 * qJ(2))) * MDP(6) + (-g(1) * (t34 * t45 - t38 * t42) - g(2) * (t34 * t42 + t38 * t45)) * MDP(17) + (-g(1) * t31 - g(2) * t29) * MDP(23) + (g(1) * t30 - g(2) * t28) * MDP(24) + (MDP(2) - MDP(4) + MDP(16)) * t32 + (-MDP(12) * t41 - MDP(14) * t35 + MDP(3) - MDP(5) - t54) * t33; (-MDP(17) - MDP(6)) * t32; (MDP(17) * pkin(3) + MDP(12)) * (g(3) * t41 - t32 * t44) + (t55 * t35 + t54) * g(3) + (t41 * MDP(13) + MDP(15) * t35 - t55 * t36) * t32; -t33 * MDP(17); (-g(1) * t28 - g(2) * t30 + t40 * t53) * MDP(23) + (g(1) * t29 - g(2) * t31 + t43 * t53) * MDP(24);];
taug = t1;
