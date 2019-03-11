% Calculate Gravitation load on the joints for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPRPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:40:03
% EndTime: 2019-03-09 01:40:04
% DurationCPUTime: 0.31s
% Computational Cost: add. (245->66), mult. (208->93), div. (0->0), fcn. (180->12), ass. (0->37)
t64 = qJ(1) + pkin(9);
t57 = sin(t64);
t60 = cos(t64);
t78 = g(1) * t60 + g(2) * t57;
t92 = MDP(15) - MDP(18);
t63 = pkin(10) + qJ(4);
t56 = sin(t63);
t59 = cos(t63);
t46 = -g(3) * t59 + t78 * t56;
t89 = g(3) * t56;
t70 = sin(qJ(1));
t87 = t70 * pkin(1);
t86 = t57 * t59;
t65 = sin(pkin(11));
t85 = t57 * t65;
t67 = cos(pkin(11));
t84 = t57 * t67;
t62 = pkin(11) + qJ(6);
t55 = sin(t62);
t83 = t60 * t55;
t58 = cos(t62);
t82 = t60 * t58;
t81 = t60 * t65;
t80 = t60 * t67;
t79 = MDP(19) + MDP(8);
t77 = g(1) * t57 - g(2) * t60;
t75 = t59 * pkin(4) + t56 * qJ(5);
t68 = cos(pkin(10));
t73 = t68 * pkin(3) + pkin(2) + t75;
t71 = cos(qJ(1));
t69 = -pkin(7) - qJ(3);
t61 = t71 * pkin(1);
t51 = t57 * t55 + t59 * t82;
t50 = t57 * t58 - t59 * t83;
t49 = -t58 * t86 + t83;
t48 = t55 * t86 + t82;
t1 = [(g(1) * t71 + g(2) * t70) * MDP(3) - t78 * MDP(7) + (-g(1) * (-t57 * pkin(2) + t60 * qJ(3) - t87) - g(2) * (t60 * pkin(2) + t57 * qJ(3) + t61)) * MDP(8) + (-g(1) * (-t59 * t84 + t81) - g(2) * (t59 * t80 + t85)) * MDP(16) + (-g(1) * (t59 * t85 + t80) - g(2) * (-t59 * t81 + t84)) * MDP(17) + (g(1) * t87 - g(2) * t61 + (g(1) * t69 - g(2) * t73) * t60 + (g(1) * t73 + g(2) * t69) * t57) * MDP(19) + (-g(1) * t49 - g(2) * t51) * MDP(25) + (-g(1) * t48 - g(2) * t50) * MDP(26) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t70 - g(2) * t71) + (t59 * MDP(14) + MDP(5) * t68 - MDP(6) * sin(pkin(10)) - t92 * t56) * t77; (-MDP(4) - t79) * g(3); -t79 * t77; (-g(3) * t75 + t78 * (pkin(4) * t56 - qJ(5) * t59)) * MDP(19) + t92 * (t78 * t59 + t89) + (MDP(16) * t67 - MDP(17) * t65 + MDP(25) * t58 - MDP(26) * t55 + MDP(14)) * t46; -t46 * MDP(19); (-g(1) * t50 + g(2) * t48 + t55 * t89) * MDP(25) + (g(1) * t51 - g(2) * t49 + t58 * t89) * MDP(26);];
taug  = t1;
