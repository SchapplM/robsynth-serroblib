% Calculate Gravitation load on the joints for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S4RRRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:41
% EndTime: 2019-12-31 17:15:42
% DurationCPUTime: 0.06s
% Computational Cost: add. (75->28), mult. (100->39), div. (0->0), fcn. (75->6), ass. (0->15)
t29 = qJ(2) + qJ(3);
t25 = sin(t29);
t26 = cos(t29);
t31 = sin(qJ(1));
t33 = cos(qJ(1));
t35 = g(1) * t33 + g(2) * t31;
t34 = -g(3) * t26 + t35 * t25;
t37 = t34 * MDP(16) + (g(3) * t25 + t35 * t26) * MDP(17);
t32 = cos(qJ(2));
t36 = t32 * pkin(2) + pkin(3) * t26;
t22 = g(1) * t31 - g(2) * t33;
t30 = sin(qJ(2));
t28 = -qJ(4) - pkin(6) - pkin(5);
t20 = pkin(1) + t36;
t1 = [(-g(1) * (-t31 * t20 - t33 * t28) - g(2) * (t33 * t20 - t31 * t28)) * MDP(19) + (MDP(3) - MDP(18)) * t35 + (-t30 * MDP(10) + MDP(16) * t26 - MDP(17) * t25 + t32 * MDP(9) + MDP(2)) * t22; (-g(3) * t32 + t35 * t30) * MDP(9) + (g(3) * t30 + t35 * t32) * MDP(10) + (-g(3) * t36 - t35 * (-t30 * pkin(2) - pkin(3) * t25)) * MDP(19) + t37; t34 * MDP(19) * pkin(3) + t37; -t22 * MDP(19);];
taug = t1;
