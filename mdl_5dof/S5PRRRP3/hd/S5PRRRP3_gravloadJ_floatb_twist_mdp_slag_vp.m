% Calculate Gravitation load on the joints for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:23:27
% EndTime: 2021-01-15 16:23:28
% DurationCPUTime: 0.11s
% Computational Cost: add. (156->32), mult. (136->40), div. (0->0), fcn. (103->6), ass. (0->18)
t48 = MDP(17) + MDP(19);
t47 = MDP(18) + MDP(20);
t41 = qJ(3) + qJ(4);
t37 = cos(t41);
t43 = cos(qJ(3));
t46 = t43 * pkin(3) + pkin(4) * t37;
t36 = sin(t41);
t40 = pkin(8) + qJ(2);
t34 = sin(t40);
t35 = cos(t40);
t44 = g(1) * t35 + g(2) * t34;
t25 = -g(3) * t37 + t44 * t36;
t45 = t47 * (g(3) * t36 + t44 * t37) + t48 * t25;
t29 = g(1) * t34 - g(2) * t35;
t42 = sin(qJ(3));
t39 = -qJ(5) - pkin(7) - pkin(6);
t31 = pkin(2) + t46;
t1 = [(-MDP(1) - MDP(22)) * g(3); (-g(1) * (-t34 * t31 - t35 * t39) - g(2) * (t35 * t31 - t34 * t39)) * MDP(22) + (MDP(4) - MDP(21)) * t44 + (t43 * MDP(10) - t42 * MDP(11) - t47 * t36 + t48 * t37 + MDP(3)) * t29; (-g(3) * t43 + t44 * t42) * MDP(10) + (g(3) * t42 + t44 * t43) * MDP(11) + (-g(3) * t46 - t44 * (-t42 * pkin(3) - pkin(4) * t36)) * MDP(22) + t45; t25 * MDP(22) * pkin(4) + t45; -t29 * MDP(22);];
taug = t1;
