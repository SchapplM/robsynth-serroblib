% Calculate Gravitation load on the joints for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PPRRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:45
% EndTime: 2019-12-31 17:35:45
% DurationCPUTime: 0.05s
% Computational Cost: add. (89->19), mult. (103->32), div. (0->0), fcn. (112->8), ass. (0->16)
t43 = sin(pkin(8));
t42 = qJ(3) + qJ(4);
t41 = cos(t42);
t40 = sin(t42);
t32 = cos(pkin(8));
t25 = -t32 * t41 - t43 * t40;
t26 = t32 * t40 - t43 * t41;
t33 = sin(qJ(5));
t35 = cos(qJ(5));
t37 = -g(1) * t25 - g(2) * t26;
t39 = t37 * MDP(8) + (t35 * MDP(14) - t33 * MDP(15) + MDP(7)) * (g(1) * t26 - g(2) * t25);
t36 = cos(qJ(3));
t34 = sin(qJ(3));
t28 = t32 * t34 - t43 * t36;
t27 = -t32 * t36 - t43 * t34;
t1 = [(-MDP(1) - MDP(2)) * g(3); (-g(1) * t43 + g(2) * t32) * MDP(2); (g(1) * t28 - g(2) * t27) * MDP(4) + (-g(1) * t27 - g(2) * t28) * MDP(5) + t39; t39; (g(3) * t35 + t37 * t33) * MDP(14) + (-g(3) * t33 + t37 * t35) * MDP(15);];
taug = t1;
