% Calculate Gravitation load on the joints for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S4RRRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:37
% EndTime: 2019-12-31 17:24:37
% DurationCPUTime: 0.05s
% Computational Cost: add. (108->22), mult. (108->30), div. (0->0), fcn. (88->8), ass. (0->14)
t26 = qJ(2) + qJ(3);
t25 = qJ(4) + t26;
t21 = sin(t25);
t22 = cos(t25);
t28 = sin(qJ(1));
t30 = cos(qJ(1));
t32 = g(1) * t30 + g(2) * t28;
t34 = (-g(3) * t22 + t32 * t21) * MDP(23) + (g(3) * t21 + t32 * t22) * MDP(24);
t23 = sin(t26);
t24 = cos(t26);
t33 = (-g(3) * t24 + t32 * t23) * MDP(16) + (g(3) * t23 + t32 * t24) * MDP(17) + t34;
t29 = cos(qJ(2));
t27 = sin(qJ(2));
t1 = [t32 * MDP(3) + (-t27 * MDP(10) + MDP(16) * t24 - MDP(17) * t23 + MDP(23) * t22 - MDP(24) * t21 + t29 * MDP(9) + MDP(2)) * (g(1) * t28 - g(2) * t30); (-g(3) * t29 + t32 * t27) * MDP(9) + (g(3) * t27 + t32 * t29) * MDP(10) + t33; t33; t34;];
taug = t1;
