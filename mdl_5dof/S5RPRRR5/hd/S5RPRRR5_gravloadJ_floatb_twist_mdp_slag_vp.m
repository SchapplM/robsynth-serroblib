% Calculate Gravitation load on the joints for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:54:08
% EndTime: 2020-01-03 11:54:08
% DurationCPUTime: 0.05s
% Computational Cost: add. (137->22), mult. (100->30), div. (0->0), fcn. (76->8), ass. (0->14)
t34 = qJ(4) + qJ(5);
t32 = sin(t34);
t33 = cos(t34);
t31 = qJ(1) + pkin(9) + qJ(3);
t29 = sin(t31);
t30 = cos(t31);
t41 = g(2) * t29 - g(3) * t30;
t43 = (-g(1) * t33 + t41 * t32) * MDP(20) + (g(1) * t32 + t41 * t33) * MDP(21);
t35 = sin(qJ(4));
t37 = cos(qJ(4));
t39 = t41 * MDP(7) + (-t37 * MDP(13) + t35 * MDP(14) - t33 * MDP(20) + t32 * MDP(21) - MDP(6)) * (g(2) * t30 + g(3) * t29);
t38 = cos(qJ(1));
t36 = sin(qJ(1));
t1 = [(g(2) * t36 - g(3) * t38) * MDP(3) + t39 + (MDP(4) * pkin(1) + MDP(2)) * (-g(2) * t38 - g(3) * t36); -g(1) * MDP(4); t39; (-g(1) * t37 + t41 * t35) * MDP(13) + (g(1) * t35 + t41 * t37) * MDP(14) + t43; t43;];
taug = t1;
