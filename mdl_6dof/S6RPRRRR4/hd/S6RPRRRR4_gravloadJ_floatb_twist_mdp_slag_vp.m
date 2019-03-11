% Calculate Gravitation load on the joints for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:05:54
% EndTime: 2019-03-09 07:05:55
% DurationCPUTime: 0.16s
% Computational Cost: add. (286->46), mult. (225->64), div. (0->0), fcn. (196->12), ass. (0->27)
t54 = pkin(11) + qJ(3);
t53 = qJ(4) + t54;
t50 = qJ(5) + t53;
t46 = sin(t50);
t70 = g(3) * t46;
t57 = sin(qJ(6));
t60 = cos(qJ(1));
t68 = t57 * t60;
t58 = sin(qJ(1));
t67 = t58 * t57;
t59 = cos(qJ(6));
t66 = t58 * t59;
t65 = t59 * t60;
t47 = cos(t50);
t63 = g(1) * t60 + g(2) * t58;
t64 = (t63 * t47 + t70) * MDP(28) + (t59 * MDP(34) - t57 * MDP(35) + MDP(27)) * (-g(3) * t47 + t63 * t46);
t44 = g(1) * t58 - g(2) * t60;
t48 = sin(t53);
t49 = cos(t53);
t62 = (-g(3) * t49 + t63 * t48) * MDP(20) + (g(3) * t48 + t63 * t49) * MDP(21) + t64;
t52 = cos(t54);
t51 = sin(t54);
t43 = t47 * t65 + t67;
t42 = -t47 * t68 + t66;
t41 = -t47 * t66 + t68;
t40 = t47 * t67 + t65;
t1 = [(-g(1) * (-t58 * pkin(1) + qJ(2) * t60) - g(2) * (pkin(1) * t60 + t58 * qJ(2))) * MDP(7) + (-g(1) * t41 - g(2) * t43) * MDP(34) + (-g(1) * t40 - g(2) * t42) * MDP(35) + (MDP(3) - MDP(6)) * t63 + (t52 * MDP(13) - t51 * MDP(14) + MDP(20) * t49 - MDP(21) * t48 + MDP(27) * t47 - MDP(28) * t46 + MDP(4) * cos(pkin(11)) - MDP(5) * sin(pkin(11)) + MDP(2)) * t44; -t44 * MDP(7); (-g(3) * t52 + t63 * t51) * MDP(13) + (g(3) * t51 + t63 * t52) * MDP(14) + t62; t62; t64; (-g(1) * t42 + g(2) * t40 + t57 * t70) * MDP(34) + (g(1) * t43 - g(2) * t41 + t59 * t70) * MDP(35);];
taug  = t1;
