% Calculate Gravitation load on the joints for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:46:22
% EndTime: 2021-01-15 10:46:23
% DurationCPUTime: 0.08s
% Computational Cost: add. (93->31), mult. (108->40), div. (0->0), fcn. (83->8), ass. (0->16)
t28 = qJ(2) + pkin(7);
t27 = qJ(4) + t28;
t22 = sin(t27);
t23 = cos(t27);
t31 = sin(qJ(1));
t33 = cos(qJ(1));
t35 = g(1) * t33 + g(2) * t31;
t36 = (-g(3) * t23 + t35 * t22) * MDP(20) + (g(3) * t22 + t35 * t23) * MDP(21);
t20 = g(1) * t31 - g(2) * t33;
t32 = cos(qJ(2));
t30 = sin(qJ(2));
t29 = -qJ(3) - pkin(5);
t26 = cos(t28);
t25 = sin(t28);
t24 = t32 * pkin(2) + pkin(1);
t1 = [(-g(1) * (-t31 * t24 - t29 * t33) - g(2) * (t33 * t24 - t31 * t29)) * MDP(14) + (MDP(3) - MDP(13)) * t35 + (-t30 * MDP(10) + t26 * MDP(11) - t25 * MDP(12) + MDP(20) * t23 - MDP(21) * t22 + t32 * MDP(9) + MDP(2)) * t20; (g(3) * t30 + t35 * t32) * MDP(10) + (-g(3) * t26 + t35 * t25) * MDP(11) + (g(3) * t25 + t35 * t26) * MDP(12) + t36 + (MDP(14) * pkin(2) + MDP(9)) * (-g(3) * t32 + t35 * t30); -t20 * MDP(14); t36;];
taug = t1;
