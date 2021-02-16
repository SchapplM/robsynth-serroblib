% Calculate Gravitation load on the joints for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:04
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:04:29
% EndTime: 2021-01-15 11:04:30
% DurationCPUTime: 0.07s
% Computational Cost: add. (111->28), mult. (120->35), div. (0->0), fcn. (89->6), ass. (0->17)
t47 = MDP(12) + MDP(14);
t46 = MDP(13) + MDP(15);
t40 = cos(qJ(3));
t33 = t40 * pkin(3) + pkin(2);
t36 = qJ(1) + qJ(2);
t34 = sin(t36);
t35 = cos(t36);
t37 = -qJ(4) - pkin(6);
t45 = t35 * t33 - t34 * t37;
t31 = g(1) * t35 + g(2) * t34;
t30 = g(1) * t34 - g(2) * t35;
t44 = -t34 * t33 - t35 * t37;
t38 = sin(qJ(3));
t43 = (-MDP(16) + MDP(6)) * t31 + (-t46 * t38 + t47 * t40 + MDP(5)) * t30;
t41 = cos(qJ(1));
t39 = sin(qJ(1));
t1 = [(g(1) * t39 - g(2) * t41) * MDP(2) + (g(1) * t41 + g(2) * t39) * MDP(3) + (-g(1) * (-t39 * pkin(1) + t44) - g(2) * (t41 * pkin(1) + t45)) * MDP(17) + t43; (-g(1) * t44 - g(2) * t45) * MDP(17) + t43; t46 * (g(3) * t38 + t31 * t40) + (MDP(17) * pkin(3) + t47) * (-g(3) * t40 + t31 * t38); -t30 * MDP(17);];
taug = t1;
