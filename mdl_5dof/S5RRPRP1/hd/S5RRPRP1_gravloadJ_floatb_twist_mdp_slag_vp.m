% Calculate Gravitation load on the joints for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:20:00
% EndTime: 2022-01-20 10:20:00
% DurationCPUTime: 0.10s
% Computational Cost: add. (170->39), mult. (138->49), div. (0->0), fcn. (99->8), ass. (0->25)
t65 = MDP(13) + MDP(15);
t64 = MDP(14) + MDP(16);
t49 = qJ(1) + qJ(2);
t46 = sin(t49);
t63 = pkin(2) * t46;
t52 = sin(qJ(1));
t62 = t52 * pkin(1);
t45 = pkin(8) + t49;
t41 = sin(t45);
t42 = cos(t45);
t47 = cos(t49);
t43 = pkin(2) * t47;
t53 = cos(qJ(4));
t44 = t53 * pkin(4) + pkin(3);
t50 = -qJ(5) - pkin(7);
t61 = -t41 * t50 + t42 * t44 + t43;
t60 = g(1) * t42 + g(2) * t41;
t59 = g(1) * t41 - g(2) * t42;
t58 = g(1) * t46 - g(2) * t47;
t51 = sin(qJ(4));
t57 = -t60 * MDP(17) + t58 * MDP(5) + (g(1) * t47 + g(2) * t46) * MDP(6) + (-t64 * t51 + t65 * t53) * t59;
t56 = -t41 * t44 - t42 * t50 - t63;
t54 = cos(qJ(1));
t48 = t54 * pkin(1);
t1 = [(g(1) * t52 - g(2) * t54) * MDP(2) + (g(1) * t54 + g(2) * t52) * MDP(3) + (-g(1) * (-t62 - t63) - g(2) * (t43 + t48)) * MDP(7) + (-g(1) * (t56 - t62) - g(2) * (t48 + t61)) * MDP(18) + t57; t58 * pkin(2) * MDP(7) + (-g(1) * t56 - g(2) * t61) * MDP(18) + t57; (-MDP(18) - MDP(7)) * g(3); t64 * (g(3) * t51 + t60 * t53) + (MDP(18) * pkin(4) + t65) * (-g(3) * t53 + t60 * t51); -t59 * MDP(18);];
taug = t1;
