% Calculate Gravitation load on the joints for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:27
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:27:06
% EndTime: 2021-01-14 22:27:06
% DurationCPUTime: 0.04s
% Computational Cost: add. (67->21), mult. (76->25), div. (0->0), fcn. (55->4), ass. (0->12)
t30 = MDP(10) + MDP(12);
t29 = MDP(11) + MDP(13);
t23 = pkin(6) + qJ(2);
t21 = sin(t23);
t22 = cos(t23);
t28 = g(1) * t22 + g(2) * t21;
t18 = g(1) * t21 - g(2) * t22;
t26 = cos(qJ(3));
t25 = sin(qJ(3));
t24 = -qJ(4) - pkin(5);
t20 = t26 * pkin(3) + pkin(2);
t1 = [(-MDP(1) - MDP(15)) * g(3); (-g(1) * (-t21 * t20 - t22 * t24) - g(2) * (t22 * t20 - t21 * t24)) * MDP(15) + (MDP(4) - MDP(14)) * t28 + (-t29 * t25 + t30 * t26 + MDP(3)) * t18; t29 * (g(3) * t25 + t28 * t26) + (MDP(15) * pkin(3) + t30) * (-g(3) * t26 + t28 * t25); -t18 * MDP(15);];
taug = t1;
