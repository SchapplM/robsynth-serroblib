% Calculate Gravitation load on the joints for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:49
% EndTime: 2019-12-31 18:17:50
% DurationCPUTime: 0.06s
% Computational Cost: add. (125->29), mult. (90->37), div. (0->0), fcn. (62->8), ass. (0->13)
t35 = qJ(1) + pkin(8);
t34 = qJ(3) + t35;
t32 = sin(t34);
t33 = cos(t34);
t43 = t33 * pkin(3) + t32 * qJ(4);
t42 = -t32 * pkin(3) + t33 * qJ(4);
t26 = g(1) * t32 - g(2) * t33;
t36 = sin(qJ(5));
t38 = cos(qJ(5));
t40 = (MDP(6) - MDP(8)) * t26 + (-t36 * MDP(16) - t38 * MDP(17) + MDP(7) - MDP(9)) * (g(1) * t33 + g(2) * t32);
t39 = cos(qJ(1));
t37 = sin(qJ(1));
t1 = [(g(1) * t39 + g(2) * t37) * MDP(3) + (-g(1) * (-pkin(2) * sin(t35) - t37 * pkin(1) + t42) - g(2) * (pkin(2) * cos(t35) + t39 * pkin(1) + t43)) * MDP(10) + t40 + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t37 - g(2) * t39); (-MDP(10) - MDP(4)) * g(3); (-g(1) * t42 - g(2) * t43) * MDP(10) + t40; -t26 * MDP(10); (g(3) * t36 - t26 * t38) * MDP(16) + (g(3) * t38 + t26 * t36) * MDP(17);];
taug = t1;
