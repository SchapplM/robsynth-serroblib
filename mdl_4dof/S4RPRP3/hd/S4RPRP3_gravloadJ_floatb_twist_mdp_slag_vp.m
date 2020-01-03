% Calculate Gravitation load on the joints for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S4RPRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:48
% EndTime: 2019-12-31 16:42:48
% DurationCPUTime: 0.10s
% Computational Cost: add. (49->24), mult. (61->34), div. (0->0), fcn. (41->6), ass. (0->12)
t17 = qJ(1) + pkin(6);
t15 = sin(t17);
t16 = cos(t17);
t26 = g(1) * t16 + g(2) * t15;
t25 = g(1) * t15 - g(2) * t16;
t22 = cos(qJ(1));
t21 = cos(qJ(3));
t20 = sin(qJ(1));
t19 = sin(qJ(3));
t18 = -qJ(4) - pkin(5);
t14 = t21 * pkin(3) + pkin(2);
t1 = [(g(1) * t22 + g(2) * t20) * MDP(3) - t26 * MDP(12) + (-g(1) * (-t20 * pkin(1) - t15 * t14 - t16 * t18) - g(2) * (t22 * pkin(1) + t16 * t14 - t15 * t18)) * MDP(13) + (MDP(10) * t21 - MDP(11) * t19) * t25 + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t20 - g(2) * t22); (-MDP(13) - MDP(4)) * g(3); (g(3) * t19 + t26 * t21) * MDP(11) + (MDP(13) * pkin(3) + MDP(10)) * (-g(3) * t21 + t26 * t19); -t25 * MDP(13);];
taug = t1;
