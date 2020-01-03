% Calculate Gravitation load on the joints for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:56
% EndTime: 2019-12-31 16:37:57
% DurationCPUTime: 0.07s
% Computational Cost: add. (58->24), mult. (60->34), div. (0->0), fcn. (42->8), ass. (0->11)
t19 = qJ(1) + pkin(6);
t15 = sin(t19);
t17 = cos(t19);
t26 = g(1) * t17 + g(2) * t15;
t25 = g(1) * t15 - g(2) * t17;
t23 = cos(qJ(1));
t22 = sin(qJ(1));
t18 = pkin(7) + qJ(4);
t16 = cos(t18);
t14 = sin(t18);
t1 = [(g(1) * t23 + g(2) * t22) * MDP(3) - t26 * MDP(7) + (-g(1) * (-t22 * pkin(1) - t15 * pkin(2) + t17 * qJ(3)) - g(2) * (t23 * pkin(1) + t17 * pkin(2) + t15 * qJ(3))) * MDP(8) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t22 - g(2) * t23) + (t16 * MDP(14) - t14 * MDP(15) + MDP(5) * cos(pkin(7)) - MDP(6) * sin(pkin(7))) * t25; (-MDP(4) - MDP(8)) * g(3); -t25 * MDP(8); (-g(3) * t16 + t26 * t14) * MDP(14) + (g(3) * t14 + t26 * t16) * MDP(15);];
taug = t1;
