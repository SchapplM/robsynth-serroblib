% Calculate Gravitation load on the joints for
% S4RRPR4
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
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RRPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:33
% EndTime: 2019-12-31 17:02:34
% DurationCPUTime: 0.08s
% Computational Cost: add. (102->26), mult. (95->35), div. (0->0), fcn. (70->8), ass. (0->14)
t39 = qJ(1) + qJ(2);
t36 = sin(t39);
t37 = cos(t39);
t46 = t37 * pkin(2) + t36 * qJ(3);
t45 = -t36 * pkin(2) + t37 * qJ(3);
t30 = g(1) * t37 + g(2) * t36;
t29 = g(1) * t36 - g(2) * t37;
t38 = pkin(7) + qJ(4);
t34 = sin(t38);
t35 = cos(t38);
t44 = (MDP(6) - MDP(9)) * t30 + (t35 * MDP(16) - t34 * MDP(17) + MDP(7) * cos(pkin(7)) - MDP(8) * sin(pkin(7)) + MDP(5)) * t29;
t43 = cos(qJ(1));
t42 = sin(qJ(1));
t1 = [(g(1) * t42 - g(2) * t43) * MDP(2) + (g(1) * t43 + g(2) * t42) * MDP(3) + (-g(1) * (-t42 * pkin(1) + t45) - g(2) * (t43 * pkin(1) + t46)) * MDP(10) + t44; (-g(1) * t45 - g(2) * t46) * MDP(10) + t44; -t29 * MDP(10); (-g(3) * t35 + t30 * t34) * MDP(16) + (g(3) * t34 + t30 * t35) * MDP(17);];
taug = t1;
