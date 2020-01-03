% Calculate Gravitation load on the joints for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPPRR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:53
% EndTime: 2019-12-31 17:59:54
% DurationCPUTime: 0.11s
% Computational Cost: add. (87->37), mult. (111->58), div. (0->0), fcn. (96->8), ass. (0->18)
t36 = cos(qJ(4));
t44 = g(3) * t36;
t32 = sin(qJ(5));
t33 = sin(qJ(4));
t43 = t32 * t33;
t35 = cos(qJ(5));
t42 = t33 * t35;
t31 = qJ(1) + pkin(8);
t29 = sin(t31);
t30 = cos(t31);
t40 = g(1) * t29 - g(2) * t30;
t37 = cos(qJ(1));
t34 = sin(qJ(1));
t27 = -t29 * t32 + t30 * t42;
t26 = t29 * t35 + t30 * t43;
t25 = t29 * t42 + t30 * t32;
t24 = -t29 * t43 + t30 * t35;
t1 = [(g(1) * t37 + g(2) * t34) * MDP(3) - t40 * MDP(5) + (-g(1) * (-t34 * pkin(1) - t29 * pkin(2) + t30 * qJ(3)) - g(2) * (t37 * pkin(1) + t30 * pkin(2) + t29 * qJ(3))) * MDP(7) + (-g(1) * t27 - g(2) * t25) * MDP(20) + (g(1) * t26 - g(2) * t24) * MDP(21) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t34 - g(2) * t37) + (t33 * MDP(13) + t36 * MDP(14) + MDP(6)) * (-g(1) * t30 - g(2) * t29); (-MDP(4) - MDP(7)) * g(3); -t40 * MDP(7); (t40 * t33 + t44) * MDP(14) + (-MDP(20) * t35 + MDP(21) * t32 - MDP(13)) * (-g(3) * t33 + t40 * t36); (-g(1) * t24 - g(2) * t26 + t32 * t44) * MDP(20) + (g(1) * t25 - g(2) * t27 + t35 * t44) * MDP(21);];
taug = t1;
