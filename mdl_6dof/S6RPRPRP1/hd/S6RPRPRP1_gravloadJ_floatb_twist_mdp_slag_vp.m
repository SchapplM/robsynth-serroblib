% Calculate Gravitation load on the joints for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S6RPRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:03:04
% EndTime: 2019-03-09 03:03:05
% DurationCPUTime: 0.32s
% Computational Cost: add. (207->63), mult. (203->86), div. (0->0), fcn. (169->10), ass. (0->35)
t63 = qJ(1) + pkin(9);
t57 = sin(t63);
t59 = cos(t63);
t78 = g(1) * t59 + g(2) * t57;
t62 = qJ(3) + pkin(10);
t56 = sin(t62);
t58 = cos(t62);
t94 = -g(3) * t58 + t78 * t56;
t89 = g(3) * t56;
t68 = sin(qJ(1));
t87 = t68 * pkin(1);
t66 = sin(qJ(5));
t86 = t57 * t66;
t69 = cos(qJ(5));
t85 = t57 * t69;
t84 = t59 * t66;
t83 = t59 * t69;
t70 = cos(qJ(3));
t60 = t70 * pkin(3);
t55 = t60 + pkin(2);
t71 = cos(qJ(1));
t82 = t71 * pkin(1) + t59 * t55;
t81 = MDP(13) + MDP(22);
t65 = -qJ(4) - pkin(7);
t79 = pkin(5) * t66 - t65;
t77 = g(1) * t57 - g(2) * t59;
t54 = t69 * pkin(5) + pkin(4);
t64 = -qJ(6) - pkin(8);
t75 = t58 * t54 - t56 * t64;
t50 = -t58 * t84 + t85;
t48 = t58 * t86 + t83;
t67 = sin(qJ(3));
t51 = t58 * t83 + t86;
t49 = -t58 * t85 + t84;
t1 = [(g(1) * t71 + g(2) * t68) * MDP(3) - t78 * MDP(12) + (-g(1) * (-t57 * t55 - t59 * t65 - t87) - g(2) * (-t57 * t65 + t82)) * MDP(13) + (-g(1) * t49 - g(2) * t51) * MDP(19) + (-g(1) * t48 - g(2) * t50) * MDP(20) + (g(1) * t87 - g(2) * t82 + (-g(1) * t79 - g(2) * t75) * t59 + (-g(1) * (-t55 - t75) - g(2) * t79) * t57) * MDP(22) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t68 - g(2) * t71) + (t70 * MDP(10) - t67 * MDP(11) + t56 * MDP(21)) * t77; (-MDP(4) - t81) * g(3); (g(3) * t67 + t78 * t70) * MDP(11) + (-t78 * t58 - t89) * MDP(21) + (-g(3) * (t60 + t75) + t78 * (pkin(3) * t67 + t54 * t56 + t58 * t64)) * MDP(22) + (pkin(3) * MDP(13) + MDP(10)) * (-g(3) * t70 + t78 * t67) + (MDP(19) * t69 - MDP(20) * t66) * t94; -t81 * t77; (g(1) * t51 - g(2) * t49 + t69 * t89) * MDP(20) + (pkin(5) * MDP(22) + MDP(19)) * (-g(1) * t50 + g(2) * t48 + t66 * t89); -t94 * MDP(22);];
taug  = t1;
