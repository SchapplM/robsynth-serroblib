% Calculate Gravitation load on the joints for
% S4PRRR3
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
%   see S4PRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4PRRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:39
% EndTime: 2019-12-31 16:31:39
% DurationCPUTime: 0.04s
% Computational Cost: add. (68->15), mult. (51->21), div. (0->0), fcn. (38->6), ass. (0->11)
t22 = pkin(7) + qJ(2);
t21 = qJ(3) + t22;
t17 = sin(t21);
t18 = cos(t21);
t23 = sin(qJ(4));
t24 = cos(qJ(4));
t26 = g(1) * t18 + g(2) * t17;
t27 = t26 * MDP(7) + (t24 * MDP(13) - t23 * MDP(14) + MDP(6)) * (g(1) * t17 - g(2) * t18);
t20 = cos(t22);
t19 = sin(t22);
t1 = [-g(3) * MDP(1); (g(1) * t19 - g(2) * t20) * MDP(3) + (g(1) * t20 + g(2) * t19) * MDP(4) + t27; t27; (-g(3) * t24 + t26 * t23) * MDP(13) + (g(3) * t23 + t26 * t24) * MDP(14);];
taug = t1;
