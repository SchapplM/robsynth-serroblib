% Calculate Gravitation load on the joints for
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PPRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4PPRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:52
% EndTime: 2019-12-31 16:19:52
% DurationCPUTime: 0.05s
% Computational Cost: add. (22->17), mult. (57->30), div. (0->0), fcn. (52->6), ass. (0->11)
t19 = cos(qJ(3));
t24 = g(3) * t19;
t16 = sin(qJ(4));
t17 = sin(qJ(3));
t23 = t16 * t17;
t18 = cos(qJ(4));
t22 = t17 * t18;
t14 = sin(pkin(6));
t15 = cos(pkin(6));
t21 = g(1) * t14 - g(2) * t15;
t1 = [(-MDP(1) - MDP(2)) * g(3); -t21 * MDP(2); (t21 * t17 + t24) * MDP(5) + (-MDP(11) * t18 + MDP(12) * t16 - MDP(4)) * (-g(3) * t17 + t21 * t19); (-g(1) * (-t14 * t23 + t15 * t18) - g(2) * (t14 * t18 + t15 * t23) + t16 * t24) * MDP(11) + (-g(1) * (-t14 * t22 - t15 * t16) - g(2) * (-t14 * t16 + t15 * t22) + t18 * t24) * MDP(12);];
taug = t1;
