% Calculate Gravitation load on the joints for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:23:36
% EndTime: 2021-01-15 21:23:37
% DurationCPUTime: 0.14s
% Computational Cost: add. (180->38), mult. (154->48), div. (0->0), fcn. (121->10), ass. (0->20)
t38 = qJ(2) + pkin(9);
t37 = qJ(4) + t38;
t33 = qJ(5) + t37;
t29 = sin(t33);
t30 = cos(t33);
t41 = sin(qJ(1));
t43 = cos(qJ(1));
t45 = g(1) * t43 + g(2) * t41;
t47 = (-g(3) * t30 + t29 * t45) * MDP(27) + (g(3) * t29 + t30 * t45) * MDP(28);
t31 = sin(t37);
t32 = cos(t37);
t46 = (-g(3) * t32 + t31 * t45) * MDP(20) + (g(3) * t31 + t32 * t45) * MDP(21) + t47;
t27 = g(1) * t41 - g(2) * t43;
t42 = cos(qJ(2));
t40 = sin(qJ(2));
t39 = -qJ(3) - pkin(6);
t36 = cos(t38);
t35 = sin(t38);
t34 = pkin(2) * t42 + pkin(1);
t1 = [(-g(1) * (-t34 * t41 - t39 * t43) - g(2) * (t34 * t43 - t39 * t41)) * MDP(14) + (MDP(3) - MDP(13)) * t45 + (-t40 * MDP(10) + t36 * MDP(11) - t35 * MDP(12) + MDP(20) * t32 - MDP(21) * t31 + MDP(27) * t30 - MDP(28) * t29 + MDP(9) * t42 + MDP(2)) * t27; (g(3) * t40 + t42 * t45) * MDP(10) + (-g(3) * t36 + t35 * t45) * MDP(11) + (g(3) * t35 + t36 * t45) * MDP(12) + t46 + (MDP(14) * pkin(2) + MDP(9)) * (-g(3) * t42 + t40 * t45); -t27 * MDP(14); t46; t47;];
taug = t1;
