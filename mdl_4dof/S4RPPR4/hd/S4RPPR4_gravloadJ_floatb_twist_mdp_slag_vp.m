% Calculate Gravitation load on the joints for
% S4RPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4RPPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:54
% EndTime: 2019-12-31 16:38:55
% DurationCPUTime: 0.07s
% Computational Cost: add. (44->22), mult. (53->32), div. (0->0), fcn. (36->6), ass. (0->9)
t15 = qJ(1) + pkin(6);
t13 = sin(t15);
t14 = cos(t15);
t21 = g(1) * t13 - g(2) * t14;
t19 = cos(qJ(1));
t18 = cos(qJ(4));
t17 = sin(qJ(1));
t16 = sin(qJ(4));
t1 = [(g(1) * t19 + g(2) * t17) * MDP(3) - t21 * MDP(5) + (-g(1) * (-t17 * pkin(1) - t13 * pkin(2) + t14 * qJ(3)) - g(2) * (t19 * pkin(1) + t14 * pkin(2) + t13 * qJ(3))) * MDP(7) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t17 - g(2) * t19) + (t16 * MDP(13) + t18 * MDP(14) + MDP(6)) * (-g(1) * t14 - g(2) * t13); (-MDP(4) - MDP(7)) * g(3); -t21 * MDP(7); (g(3) * t16 - t21 * t18) * MDP(13) + (g(3) * t18 + t21 * t16) * MDP(14);];
taug = t1;
