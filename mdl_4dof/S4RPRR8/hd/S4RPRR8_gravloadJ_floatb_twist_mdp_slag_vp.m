% Calculate Gravitation load on the joints for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S4RPRR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:15
% EndTime: 2019-12-31 16:55:15
% DurationCPUTime: 0.06s
% Computational Cost: add. (52->22), mult. (78->30), div. (0->0), fcn. (60->6), ass. (0->10)
t22 = sin(qJ(1));
t24 = cos(qJ(1));
t15 = g(1) * t22 - g(2) * t24;
t20 = qJ(3) + qJ(4);
t17 = sin(t20);
t18 = cos(t20);
t26 = (g(3) * t17 - t15 * t18) * MDP(19) + (g(3) * t18 + t15 * t17) * MDP(20);
t23 = cos(qJ(3));
t21 = sin(qJ(3));
t1 = [(-g(1) * (-t22 * pkin(1) + t24 * qJ(2)) - g(2) * (t24 * pkin(1) + t22 * qJ(2))) * MDP(6) + (MDP(2) - MDP(4)) * t15 + (-t21 * MDP(12) - t23 * MDP(13) - MDP(19) * t17 - MDP(20) * t18 + MDP(3) - MDP(5)) * (g(1) * t24 + g(2) * t22); -t15 * MDP(6); (g(3) * t21 - t15 * t23) * MDP(12) + (g(3) * t23 + t15 * t21) * MDP(13) + t26; t26;];
taug = t1;
