% Calculate Gravitation load on the joints for
% S4RPRR3
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
%   see S4RPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RPRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:18
% EndTime: 2019-12-31 16:49:19
% DurationCPUTime: 0.06s
% Computational Cost: add. (66->20), mult. (68->29), div. (0->0), fcn. (52->8), ass. (0->13)
t23 = qJ(3) + qJ(4);
t20 = sin(t23);
t21 = cos(t23);
t22 = qJ(1) + pkin(7);
t18 = sin(t22);
t19 = cos(t22);
t30 = g(1) * t19 + g(2) * t18;
t31 = (-g(3) * t21 + t30 * t20) * MDP(17) + (g(3) * t20 + t30 * t21) * MDP(18);
t27 = cos(qJ(1));
t26 = cos(qJ(3));
t25 = sin(qJ(1));
t24 = sin(qJ(3));
t1 = [(g(1) * t27 + g(2) * t25) * MDP(3) + (MDP(4) * pkin(1) + MDP(2)) * (g(1) * t25 - g(2) * t27) + (t26 * MDP(10) - t24 * MDP(11) + MDP(17) * t21 - MDP(18) * t20) * (g(1) * t18 - g(2) * t19); -g(3) * MDP(4); (-g(3) * t26 + t30 * t24) * MDP(10) + (g(3) * t24 + t30 * t26) * MDP(11) + t31; t31;];
taug = t1;
