% Calculate Gravitation load on the joints for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4PRRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:46
% EndTime: 2019-12-31 16:33:46
% DurationCPUTime: 0.09s
% Computational Cost: add. (68->21), mult. (95->33), div. (0->0), fcn. (86->8), ass. (0->17)
t27 = qJ(2) + qJ(3);
t25 = sin(t27);
t42 = g(3) * t25;
t28 = sin(pkin(7));
t30 = sin(qJ(4));
t40 = t28 * t30;
t32 = cos(qJ(4));
t39 = t28 * t32;
t29 = cos(pkin(7));
t38 = t29 * t30;
t37 = t29 * t32;
t26 = cos(t27);
t35 = g(1) * t29 + g(2) * t28;
t36 = (t35 * t26 + t42) * MDP(7) + (t32 * MDP(13) - t30 * MDP(14) + MDP(6)) * (-g(3) * t26 + t35 * t25);
t33 = cos(qJ(2));
t31 = sin(qJ(2));
t1 = [-g(3) * MDP(1); (-g(3) * t33 + t35 * t31) * MDP(3) + (g(3) * t31 + t35 * t33) * MDP(4) + t36; t36; (-g(1) * (-t26 * t38 + t39) - g(2) * (-t26 * t40 - t37) + t30 * t42) * MDP(13) + (-g(1) * (-t26 * t37 - t40) - g(2) * (-t26 * t39 + t38) + t32 * t42) * MDP(14);];
taug = t1;
