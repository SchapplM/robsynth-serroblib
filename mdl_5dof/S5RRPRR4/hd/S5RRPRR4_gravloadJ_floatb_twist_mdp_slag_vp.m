% Calculate Gravitation load on the joints for
% S5RRPRR4
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
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:33
% EndTime: 2022-01-20 10:48:33
% DurationCPUTime: 0.06s
% Computational Cost: add. (137->29), mult. (107->43), div. (0->0), fcn. (80->10), ass. (0->18)
t39 = qJ(4) + qJ(5);
t35 = sin(t39);
t37 = cos(t39);
t40 = qJ(1) + qJ(2);
t34 = pkin(9) + t40;
t32 = sin(t34);
t33 = cos(t34);
t48 = g(1) * t33 + g(2) * t32;
t49 = (-g(3) * t37 + t48 * t35) * MDP(20) + (g(3) * t35 + t48 * t37) * MDP(21);
t36 = sin(t40);
t38 = cos(t40);
t46 = g(1) * t36 - g(2) * t38;
t41 = sin(qJ(4));
t43 = cos(qJ(4));
t45 = t46 * MDP(5) + (g(1) * t38 + g(2) * t36) * MDP(6) + (t43 * MDP(13) - t41 * MDP(14) + t37 * MDP(20) - t35 * MDP(21)) * (g(1) * t32 - g(2) * t33);
t44 = cos(qJ(1));
t42 = sin(qJ(1));
t1 = [(g(1) * t42 - g(2) * t44) * MDP(2) + (g(1) * t44 + g(2) * t42) * MDP(3) + (-g(1) * (-t42 * pkin(1) - pkin(2) * t36) - g(2) * (t44 * pkin(1) + pkin(2) * t38)) * MDP(7) + t45; t46 * MDP(7) * pkin(2) + t45; -g(3) * MDP(7); (-g(3) * t43 + t48 * t41) * MDP(13) + (g(3) * t41 + t48 * t43) * MDP(14) + t49; t49;];
taug = t1;
