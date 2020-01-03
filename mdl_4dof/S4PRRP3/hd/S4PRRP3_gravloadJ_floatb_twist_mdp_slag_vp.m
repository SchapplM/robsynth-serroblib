% Calculate Gravitation load on the joints for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:56
% EndTime: 2019-12-31 16:26:56
% DurationCPUTime: 0.05s
% Computational Cost: add. (49->19), mult. (54->25), div. (0->0), fcn. (37->4), ass. (0->10)
t16 = pkin(6) + qJ(2);
t14 = sin(t16);
t15 = cos(t16);
t21 = g(1) * t15 + g(2) * t14;
t11 = g(1) * t14 - g(2) * t15;
t19 = cos(qJ(3));
t18 = sin(qJ(3));
t17 = -qJ(4) - pkin(5);
t13 = t19 * pkin(3) + pkin(2);
t1 = [(-MDP(1) - MDP(13)) * g(3); (-g(1) * (-t14 * t13 - t15 * t17) - g(2) * (t15 * t13 - t14 * t17)) * MDP(13) + (MDP(4) - MDP(12)) * t21 + (t19 * MDP(10) - t18 * MDP(11) + MDP(3)) * t11; (g(3) * t18 + t21 * t19) * MDP(11) + (MDP(13) * pkin(3) + MDP(10)) * (-g(3) * t19 + t21 * t18); -t11 * MDP(13);];
taug = t1;
