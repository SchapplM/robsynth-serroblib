% Calculate Gravitation load on the joints for
% S4PRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4PRPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:18
% EndTime: 2019-12-31 16:23:18
% DurationCPUTime: 0.07s
% Computational Cost: add. (37->20), mult. (66->35), div. (0->0), fcn. (57->8), ass. (0->16)
t21 = qJ(2) + pkin(7);
t19 = sin(t21);
t35 = g(3) * t19;
t22 = sin(pkin(6));
t24 = sin(qJ(4));
t33 = t22 * t24;
t26 = cos(qJ(4));
t32 = t22 * t26;
t23 = cos(pkin(6));
t31 = t23 * t24;
t30 = t23 * t26;
t29 = g(1) * t23 + g(2) * t22;
t27 = cos(qJ(2));
t25 = sin(qJ(2));
t20 = cos(t21);
t1 = [(-MDP(1) - MDP(5)) * g(3); (g(3) * t25 + t29 * t27) * MDP(4) + (MDP(5) * pkin(2) + MDP(3)) * (-g(3) * t27 + t29 * t25) + (MDP(11) * t26 - MDP(12) * t24) * (-g(3) * t20 + t29 * t19); (-g(1) * t22 + g(2) * t23) * MDP(5); (-g(1) * (-t20 * t31 + t32) - g(2) * (-t20 * t33 - t30) + t24 * t35) * MDP(11) + (-g(1) * (-t20 * t30 - t33) - g(2) * (-t20 * t32 + t31) + t26 * t35) * MDP(12);];
taug = t1;
