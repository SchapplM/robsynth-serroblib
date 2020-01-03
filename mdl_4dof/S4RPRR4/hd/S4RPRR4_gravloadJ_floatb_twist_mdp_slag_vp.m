% Calculate Gravitation load on the joints for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RPRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:36
% EndTime: 2019-12-31 16:50:36
% DurationCPUTime: 0.07s
% Computational Cost: add. (66->28), mult. (92->47), div. (0->0), fcn. (84->8), ass. (0->18)
t30 = sin(qJ(3));
t42 = g(3) * t30;
t29 = sin(qJ(4));
t33 = cos(qJ(3));
t40 = t29 * t33;
t32 = cos(qJ(4));
t39 = t32 * t33;
t28 = qJ(1) + pkin(7);
t26 = sin(t28);
t27 = cos(t28);
t38 = g(1) * t27 + g(2) * t26;
t34 = cos(qJ(1));
t31 = sin(qJ(1));
t25 = t26 * t29 + t27 * t39;
t24 = t26 * t32 - t27 * t40;
t23 = -t26 * t39 + t27 * t29;
t22 = t26 * t40 + t27 * t32;
t1 = [(g(1) * t34 + g(2) * t31) * MDP(3) + (-g(1) * t23 - g(2) * t25) * MDP(17) + (-g(1) * t22 - g(2) * t24) * MDP(18) + (t33 * MDP(10) - t30 * MDP(11)) * (g(1) * t26 - g(2) * t27) + (MDP(4) * pkin(1) + MDP(2)) * (g(1) * t31 - g(2) * t34); -g(3) * MDP(4); (t38 * t33 + t42) * MDP(11) + (MDP(17) * t32 - MDP(18) * t29 + MDP(10)) * (-g(3) * t33 + t38 * t30); (-g(1) * t24 + g(2) * t22 + t29 * t42) * MDP(17) + (g(1) * t25 - g(2) * t23 + t32 * t42) * MDP(18);];
taug = t1;
