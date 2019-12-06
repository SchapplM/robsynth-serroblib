% Calculate Gravitation load on the joints for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:44
% EndTime: 2019-12-05 18:01:44
% DurationCPUTime: 0.07s
% Computational Cost: add. (128->31), mult. (97->39), div. (0->0), fcn. (67->8), ass. (0->16)
t33 = qJ(1) + pkin(8);
t32 = qJ(3) + t33;
t29 = sin(t32);
t30 = cos(t32);
t37 = cos(qJ(4));
t31 = t37 * pkin(4) + pkin(3);
t34 = -qJ(5) - pkin(7);
t43 = t29 * t34 - t30 * t31;
t26 = g(2) * t29 - g(3) * t30;
t27 = g(2) * t30 + g(3) * t29;
t35 = sin(qJ(4));
t42 = (MDP(15) - MDP(7)) * t26 + (t37 * MDP(13) - t35 * MDP(14) + MDP(6)) * t27;
t40 = -t29 * t31 - t30 * t34;
t38 = cos(qJ(1));
t36 = sin(qJ(1));
t1 = [(-g(2) * t36 + g(3) * t38) * MDP(3) + (-g(2) * (-pkin(2) * cos(t33) - t38 * pkin(1) + t43) - g(3) * (-pkin(2) * sin(t33) - t36 * pkin(1) + t40)) * MDP(16) + t42 + (pkin(1) * MDP(4) + MDP(2)) * (g(2) * t38 + g(3) * t36); (-MDP(16) - MDP(4)) * g(1); (-g(2) * t43 - g(3) * t40) * MDP(16) + t42; (g(1) * t35 - t26 * t37) * MDP(14) + (MDP(16) * pkin(4) + MDP(13)) * (-g(1) * t37 - t26 * t35); -t27 * MDP(16);];
taug = t1;
