% Calculate Gravitation load on the joints for
% S5PRRPR3
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
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:42:10
% EndTime: 2021-01-15 15:42:10
% DurationCPUTime: 0.09s
% Computational Cost: add. (136->33), mult. (110->41), div. (0->0), fcn. (83->8), ass. (0->17)
t32 = qJ(3) + pkin(9);
t30 = qJ(5) + t32;
t23 = sin(t30);
t24 = cos(t30);
t31 = pkin(8) + qJ(2);
t26 = sin(t31);
t28 = cos(t31);
t37 = g(1) * t28 + g(2) * t26;
t38 = (-g(3) * t24 + t37 * t23) * MDP(21) + (g(3) * t23 + t37 * t24) * MDP(22);
t21 = g(1) * t26 - g(2) * t28;
t35 = cos(qJ(3));
t34 = sin(qJ(3));
t33 = -qJ(4) - pkin(6);
t29 = cos(t32);
t27 = sin(t32);
t25 = t35 * pkin(3) + pkin(2);
t1 = [(-MDP(1) - MDP(15)) * g(3); (-g(1) * (-t26 * t25 - t28 * t33) - g(2) * (t28 * t25 - t26 * t33)) * MDP(15) + (MDP(4) - MDP(14)) * t37 + (t35 * MDP(10) - t34 * MDP(11) + t29 * MDP(12) - t27 * MDP(13) + MDP(21) * t24 - MDP(22) * t23 + MDP(3)) * t21; (g(3) * t34 + t37 * t35) * MDP(11) + (-g(3) * t29 + t37 * t27) * MDP(12) + (g(3) * t27 + t37 * t29) * MDP(13) + t38 + (MDP(15) * pkin(3) + MDP(10)) * (-g(3) * t35 + t37 * t34); -t21 * MDP(15); t38;];
taug = t1;
