% Calculate Gravitation load on the joints for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:35:30
% EndTime: 2019-03-09 03:35:30
% DurationCPUTime: 0.15s
% Computational Cost: add. (205->46), mult. (171->66), div. (0->0), fcn. (145->10), ass. (0->27)
t46 = qJ(3) + pkin(11) + qJ(5);
t41 = sin(t46);
t66 = g(3) * t41;
t47 = qJ(1) + pkin(10);
t44 = sin(t47);
t49 = sin(qJ(6));
t64 = t44 * t49;
t52 = cos(qJ(6));
t63 = t44 * t52;
t45 = cos(t47);
t62 = t45 * t49;
t61 = t45 * t52;
t42 = cos(t46);
t59 = g(1) * t45 + g(2) * t44;
t60 = (t42 * t59 + t66) * MDP(20) + (t52 * MDP(26) - t49 * MDP(27) + MDP(19)) * (-g(3) * t42 + t41 * t59);
t58 = g(1) * t44 - g(2) * t45;
t54 = cos(qJ(1));
t53 = cos(qJ(3));
t51 = sin(qJ(1));
t50 = sin(qJ(3));
t48 = -qJ(4) - pkin(7);
t43 = t53 * pkin(3) + pkin(2);
t40 = t42 * t61 + t64;
t39 = -t42 * t62 + t63;
t38 = -t42 * t63 + t62;
t37 = t42 * t64 + t61;
t1 = [(g(1) * t54 + g(2) * t51) * MDP(3) - t59 * MDP(12) + (-g(1) * (-t51 * pkin(1) - t44 * t43 - t45 * t48) - g(2) * (t54 * pkin(1) + t45 * t43 - t44 * t48)) * MDP(13) + (-g(1) * t38 - g(2) * t40) * MDP(26) + (-g(1) * t37 - g(2) * t39) * MDP(27) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t51 - g(2) * t54) + (MDP(10) * t53 - MDP(11) * t50 + MDP(19) * t42 - MDP(20) * t41) * t58; (-MDP(13) - MDP(4)) * g(3); (g(3) * t50 + t53 * t59) * MDP(11) + t60 + (MDP(13) * pkin(3) + MDP(10)) * (-g(3) * t53 + t50 * t59); -t58 * MDP(13); t60; (-g(1) * t39 + g(2) * t37 + t49 * t66) * MDP(26) + (g(1) * t40 - g(2) * t38 + t52 * t66) * MDP(27);];
taug  = t1;
