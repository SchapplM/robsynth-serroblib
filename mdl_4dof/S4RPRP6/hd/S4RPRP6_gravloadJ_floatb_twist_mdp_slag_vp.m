% Calculate Gravitation load on the joints for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RPRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:27:41
% EndTime: 2021-01-15 10:27:42
% DurationCPUTime: 0.06s
% Computational Cost: add. (49->26), mult. (93->32), div. (0->0), fcn. (67->4), ass. (0->11)
t29 = MDP(12) + MDP(14);
t28 = MDP(13) + MDP(15);
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t19 = g(1) * t26 + g(2) * t24;
t18 = g(1) * t24 - g(2) * t26;
t25 = cos(qJ(3));
t23 = sin(qJ(3));
t22 = pkin(1) + pkin(5) + qJ(4);
t20 = t23 * pkin(3) + qJ(2);
t1 = [(-g(1) * (-t24 * pkin(1) + t26 * qJ(2)) - g(2) * (t26 * pkin(1) + t24 * qJ(2))) * MDP(6) + (-g(1) * (t20 * t26 - t22 * t24) - g(2) * (t20 * t24 + t22 * t26)) * MDP(17) + (MDP(2) - MDP(4) + MDP(16)) * t18 + (-t29 * t23 - t28 * t25 + MDP(3) - MDP(5)) * t19; (-MDP(17) - MDP(6)) * t18; t28 * (g(3) * t25 + t18 * t23) + (MDP(17) * pkin(3) + t29) * (g(3) * t23 - t18 * t25); -t19 * MDP(17);];
taug = t1;
