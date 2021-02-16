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
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:30:23
% EndTime: 2021-01-15 14:30:23
% DurationCPUTime: 0.08s
% Computational Cost: add. (107->30), mult. (134->39), div. (0->0), fcn. (103->6), ass. (0->17)
t46 = MDP(16) + MDP(18);
t45 = MDP(17) + MDP(19);
t37 = qJ(2) + qJ(3);
t34 = cos(t37);
t40 = cos(qJ(2));
t44 = t40 * pkin(2) + pkin(3) * t34;
t33 = sin(t37);
t39 = sin(qJ(1));
t41 = cos(qJ(1));
t42 = g(1) * t41 + g(2) * t39;
t24 = -g(3) * t34 + t42 * t33;
t43 = t45 * (g(3) * t33 + t42 * t34) + t46 * t24;
t30 = g(1) * t39 - g(2) * t41;
t38 = sin(qJ(2));
t36 = -qJ(4) - pkin(6) - pkin(5);
t28 = pkin(1) + t44;
t1 = [(-g(1) * (-t39 * t28 - t41 * t36) - g(2) * (t41 * t28 - t39 * t36)) * MDP(21) + (MDP(3) - MDP(20)) * t42 + (-t38 * MDP(10) + t40 * MDP(9) - t45 * t33 + t46 * t34 + MDP(2)) * t30; (-g(3) * t40 + t42 * t38) * MDP(9) + (g(3) * t38 + t42 * t40) * MDP(10) + (-g(3) * t44 - t42 * (-t38 * pkin(2) - pkin(3) * t33)) * MDP(21) + t43; t24 * MDP(21) * pkin(3) + t43; -t30 * MDP(21);];
taug = t1;
