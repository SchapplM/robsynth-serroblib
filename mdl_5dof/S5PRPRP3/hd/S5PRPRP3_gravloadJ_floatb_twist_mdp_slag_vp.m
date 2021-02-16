% Calculate Gravitation load on the joints for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:14:01
% EndTime: 2021-01-15 15:14:02
% DurationCPUTime: 0.18s
% Computational Cost: add. (108->35), mult. (160->49), div. (0->0), fcn. (139->8), ass. (0->22)
t61 = MDP(11) + MDP(13);
t60 = MDP(12) + MDP(14);
t39 = sin(pkin(7));
t40 = cos(pkin(7));
t49 = g(1) * t40 + g(2) * t39;
t38 = qJ(2) + pkin(8);
t36 = sin(t38);
t37 = cos(t38);
t59 = -g(3) * t37 + t49 * t36;
t56 = g(3) * t36;
t42 = sin(qJ(4));
t54 = t39 * t42;
t44 = cos(qJ(4));
t53 = t39 * t44;
t52 = t40 * t42;
t51 = t40 * t44;
t50 = MDP(16) + MDP(5);
t45 = cos(qJ(2));
t43 = sin(qJ(2));
t41 = -qJ(5) - pkin(6);
t35 = t44 * pkin(4) + pkin(3);
t1 = [(-MDP(1) - t50) * g(3); (g(3) * t43 + t49 * t45) * MDP(4) + (-t49 * t37 - t56) * MDP(15) + (-g(3) * (t45 * pkin(2) + t37 * t35 - t36 * t41) + t49 * (pkin(2) * t43 + t35 * t36 + t37 * t41)) * MDP(16) + (pkin(2) * MDP(5) + MDP(3)) * (-g(3) * t45 + t49 * t43) + (-t60 * t42 + t61 * t44) * t59; t50 * (-g(1) * t39 + g(2) * t40); t60 * (-g(1) * (-t37 * t51 - t54) - g(2) * (-t37 * t53 + t52) + t44 * t56) + (MDP(16) * pkin(4) + t61) * (-g(1) * (-t37 * t52 + t53) - g(2) * (-t37 * t54 - t51) + t42 * t56); -t59 * MDP(16);];
taug = t1;
