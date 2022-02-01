% Calculate Gravitation load on the joints for
% S5RRPRR5
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
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:55
% EndTime: 2022-01-20 11:02:56
% DurationCPUTime: 0.07s
% Computational Cost: add. (169->32), mult. (129->42), div. (0->0), fcn. (98->9), ass. (0->18)
t50 = qJ(1) + qJ(2);
t47 = sin(t50);
t48 = cos(t50);
t38 = g(1) * t48 + g(2) * t47;
t49 = pkin(9) + qJ(4);
t46 = qJ(5) + t49;
t41 = sin(t46);
t42 = cos(t46);
t57 = (-g(3) * t42 + t38 * t41) * MDP(22) + (g(3) * t41 + t38 * t42) * MDP(23);
t56 = t48 * pkin(2) + t47 * qJ(3);
t55 = -t47 * pkin(2) + t48 * qJ(3);
t37 = g(1) * t47 - g(2) * t48;
t44 = sin(t49);
t45 = cos(t49);
t54 = (MDP(6) - MDP(8)) * t38 + (t45 * MDP(15) - t44 * MDP(16) + t42 * MDP(22) - t41 * MDP(23) + MDP(7) * cos(pkin(9)) + MDP(5)) * t37;
t53 = cos(qJ(1));
t52 = sin(qJ(1));
t1 = [(g(1) * t52 - g(2) * t53) * MDP(2) + (g(1) * t53 + g(2) * t52) * MDP(3) + (-g(1) * (-t52 * pkin(1) + t55) - g(2) * (t53 * pkin(1) + t56)) * MDP(9) + t54; (-g(1) * t55 - g(2) * t56) * MDP(9) + t54; -t37 * MDP(9); (-g(3) * t45 + t38 * t44) * MDP(15) + (g(3) * t44 + t38 * t45) * MDP(16) + t57; t57;];
taug = t1;
