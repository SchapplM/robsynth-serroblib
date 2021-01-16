% Calculate Gravitation load on the joints for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:20
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:20:37
% EndTime: 2021-01-15 10:20:37
% DurationCPUTime: 0.10s
% Computational Cost: add. (67->26), mult. (83->34), div. (0->0), fcn. (59->6), ass. (0->14)
t35 = MDP(10) + MDP(12);
t34 = MDP(11) + MDP(13);
t24 = qJ(1) + pkin(6);
t22 = sin(t24);
t23 = cos(t24);
t33 = g(1) * t23 + g(2) * t22;
t32 = g(1) * t22 - g(2) * t23;
t29 = cos(qJ(1));
t28 = cos(qJ(3));
t27 = sin(qJ(1));
t26 = sin(qJ(3));
t25 = -qJ(4) - pkin(5);
t21 = t28 * pkin(3) + pkin(2);
t1 = [(g(1) * t29 + g(2) * t27) * MDP(3) - t33 * MDP(14) + (-g(1) * (-t27 * pkin(1) - t22 * t21 - t23 * t25) - g(2) * (t29 * pkin(1) + t23 * t21 - t22 * t25)) * MDP(15) + (-t34 * t26 + t35 * t28) * t32 + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t27 - g(2) * t29); (-MDP(15) - MDP(4)) * g(3); t34 * (g(3) * t26 + t33 * t28) + (MDP(15) * pkin(3) + t35) * (-g(3) * t28 + t33 * t26); -t32 * MDP(15);];
taug = t1;
