% Calculate Gravitation load on the joints for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:23
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:22:40
% EndTime: 2021-01-15 15:22:41
% DurationCPUTime: 0.19s
% Computational Cost: add. (166->42), mult. (147->50), div. (0->0), fcn. (111->6), ass. (0->20)
t51 = MDP(12) + MDP(16);
t50 = MDP(13) - MDP(18);
t39 = pkin(7) + qJ(2);
t34 = sin(t39);
t36 = cos(t39);
t30 = g(1) * t36 + g(2) * t34;
t47 = -MDP(15) - MDP(19);
t29 = g(1) * t34 - g(2) * t36;
t40 = qJ(3) + pkin(8);
t35 = sin(t40);
t37 = cos(t40);
t46 = t37 * pkin(4) + t35 * qJ(5);
t43 = cos(qJ(3));
t42 = sin(qJ(3));
t41 = -qJ(4) - pkin(6);
t38 = t43 * pkin(3);
t33 = t38 + pkin(2);
t31 = t36 * t33;
t25 = -g(3) * t37 + t30 * t35;
t1 = [(-MDP(1) + t47) * g(3); (-g(1) * (-t34 * t33 - t36 * t41) - g(2) * (-t34 * t41 + t31)) * MDP(15) + (-g(2) * t31 + (g(1) * t41 - g(2) * t46) * t36 + (-g(1) * (-t33 - t46) + g(2) * t41) * t34) * MDP(19) + (MDP(4) - MDP(14) - MDP(17)) * t30 + (t43 * MDP(10) - t42 * MDP(11) - t50 * t35 + t51 * t37 + MDP(3)) * t29; (g(3) * t42 + t30 * t43) * MDP(11) + (-g(3) * (t38 + t46) + t30 * (pkin(3) * t42 + pkin(4) * t35 - qJ(5) * t37)) * MDP(19) + (pkin(3) * MDP(15) + MDP(10)) * (-g(3) * t43 + t30 * t42) + t50 * (g(3) * t35 + t30 * t37) + t51 * t25; t47 * t29; -t25 * MDP(19);];
taug = t1;
