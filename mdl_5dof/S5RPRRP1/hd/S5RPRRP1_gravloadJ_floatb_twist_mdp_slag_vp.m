% Calculate Gravitation load on the joints for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:58
% EndTime: 2019-12-05 17:59:58
% DurationCPUTime: 0.12s
% Computational Cost: add. (88->37), mult. (123->47), div. (0->0), fcn. (91->6), ass. (0->16)
t35 = sin(qJ(1));
t37 = cos(qJ(1));
t24 = g(1) * t35 - g(2) * t37;
t33 = qJ(3) + qJ(4);
t26 = sin(t33);
t27 = cos(t33);
t38 = g(3) * t26 - t24 * t27;
t40 = t38 * MDP(19) + (g(3) * t27 + t24 * t26) * MDP(20);
t39 = t37 * pkin(1) + t35 * qJ(2);
t25 = g(1) * t37 + g(2) * t35;
t36 = cos(qJ(3));
t34 = sin(qJ(3));
t32 = -qJ(5) - pkin(7) - pkin(6);
t29 = t37 * qJ(2);
t22 = t34 * pkin(3) + pkin(4) * t26;
t1 = [(-g(1) * (-t35 * pkin(1) + t29) - g(2) * t39) * MDP(6) + (-g(1) * (t37 * t22 + t29 + (-pkin(1) + t32) * t35) - g(2) * (t35 * t22 - t37 * t32 + t39)) * MDP(22) + (MDP(2) - MDP(4) + MDP(21)) * t24 + (-t34 * MDP(12) - t36 * MDP(13) - MDP(19) * t26 - MDP(20) * t27 + MDP(3) - MDP(5)) * t25; (-MDP(22) - MDP(6)) * t24; (g(3) * t34 - t24 * t36) * MDP(12) + (g(3) * t36 + t24 * t34) * MDP(13) + (g(3) * t22 - t24 * (t36 * pkin(3) + pkin(4) * t27)) * MDP(22) + t40; t38 * MDP(22) * pkin(4) + t40; -t25 * MDP(22);];
taug = t1;
