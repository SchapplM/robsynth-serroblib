% Calculate Gravitation load on the joints for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:42:39
% EndTime: 2019-03-09 02:42:40
% DurationCPUTime: 0.29s
% Computational Cost: add. (199->60), mult. (192->82), div. (0->0), fcn. (159->10), ass. (0->32)
t58 = qJ(3) + pkin(10);
t52 = sin(t58);
t54 = cos(t58);
t86 = pkin(4) * t54 + qJ(5) * t52;
t59 = qJ(1) + pkin(9);
t53 = sin(t59);
t55 = cos(t59);
t74 = g(1) * t55 + g(2) * t53;
t82 = g(3) * t54;
t61 = sin(qJ(6));
t81 = t53 * t61;
t64 = cos(qJ(6));
t80 = t53 * t64;
t79 = t55 * t61;
t78 = t55 * t64;
t65 = cos(qJ(3));
t56 = t65 * pkin(3);
t51 = t56 + pkin(2);
t66 = cos(qJ(1));
t77 = t66 * pkin(1) + t55 * t51;
t75 = MDP(13) + MDP(17);
t73 = g(1) * t53 - g(2) * t55;
t60 = -qJ(4) - pkin(7);
t63 = sin(qJ(1));
t71 = -pkin(1) * t63 - t55 * t60;
t62 = sin(qJ(3));
t47 = -t52 * t81 + t78;
t46 = t52 * t80 + t79;
t45 = t52 * t79 + t80;
t44 = t52 * t78 - t81;
t43 = -t74 * t52 + t82;
t1 = [(g(1) * t66 + g(2) * t63) * MDP(3) + (-g(1) * (-t51 * t53 + t71) - g(2) * (-t53 * t60 + t77)) * MDP(13) + (-g(1) * t71 - g(2) * (t86 * t55 + t77) + (-g(1) * (-t51 - t86) + g(2) * t60) * t53) * MDP(17) + (-g(1) * t47 - g(2) * t45) * MDP(23) + (g(1) * t46 - g(2) * t44) * MDP(24) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t63 - g(2) * t66) - (MDP(12) + MDP(14)) * t74 + (t65 * MDP(10) - t62 * MDP(11) - t54 * MDP(15) + t52 * MDP(16)) * t73; (-MDP(4) - t75) * g(3); (g(3) * t62 + t65 * t74) * MDP(11) + t43 * MDP(15) + (-g(3) * (t56 + t86) + t74 * (pkin(3) * t62 + pkin(4) * t52 - qJ(5) * t54)) * MDP(17) + (pkin(3) * MDP(13) + MDP(10)) * (-g(3) * t65 + t74 * t62) + (MDP(23) * t61 + MDP(24) * t64 + MDP(16)) * (-g(3) * t52 - t74 * t54); -t75 * t73; t43 * MDP(17); (-g(1) * t44 - g(2) * t46 + t64 * t82) * MDP(23) + (g(1) * t45 - g(2) * t47 - t61 * t82) * MDP(24);];
taug  = t1;
