% Calculate Gravitation load on the joints for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPPRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPPRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:59
% EndTime: 2019-03-09 01:34:00
% DurationCPUTime: 0.28s
% Computational Cost: add. (146->52), mult. (242->74), div. (0->0), fcn. (260->10), ass. (0->29)
t58 = pkin(10) + qJ(5);
t51 = sin(t58);
t82 = MDP(20) * t51;
t61 = sin(qJ(6));
t62 = cos(qJ(6));
t81 = MDP(26) * t62 - MDP(27) * t61 + MDP(19);
t80 = g(3) * t51;
t79 = cos(qJ(1));
t78 = sin(qJ(1));
t52 = cos(t58);
t77 = t52 * t61;
t76 = t52 * t62;
t75 = t79 * pkin(1) + t78 * qJ(2);
t74 = cos(pkin(9));
t73 = sin(pkin(9));
t70 = MDP(13) + MDP(9);
t69 = t79 * pkin(2) + t75;
t68 = -t78 * pkin(1) + t79 * qJ(2);
t44 = -t78 * t73 - t79 * t74;
t45 = t79 * t73 - t78 * t74;
t67 = g(1) * t45 - g(2) * t44;
t66 = g(1) * t44 + g(2) * t45;
t65 = t44 * t61 + t45 * t76;
t64 = -t44 * t62 + t45 * t77;
t63 = -t78 * pkin(2) + t68;
t46 = g(1) * t78 - g(2) * t79;
t42 = -t44 * t76 + t45 * t61;
t41 = t44 * t77 + t45 * t62;
t1 = [(-g(1) * t68 - g(2) * t75) * MDP(6) + (-g(1) * t63 - g(2) * t69) * MDP(9) + (-g(1) * (t45 * pkin(3) + t44 * qJ(4) + t63) - g(2) * (-pkin(3) * t44 + t45 * qJ(4) + t69)) * MDP(13) + (-g(1) * t65 - g(2) * t42) * MDP(26) + (g(1) * t64 - g(2) * t41) * MDP(27) + (MDP(3) - MDP(5)) * (g(1) * t79 + g(2) * t78) + (MDP(2) + MDP(4)) * t46 + (MDP(8) - MDP(12)) * t66 + (-MDP(10) * cos(pkin(10)) + MDP(11) * sin(pkin(10)) - t52 * MDP(19) - MDP(7) + t82) * t67; (-MDP(6) - t70) * t46; t70 * g(3); -t67 * MDP(13); (t81 * t52 - t82) * g(3) + (-MDP(20) * t52 - t81 * t51) * t66; (-g(1) * t41 - g(2) * t64 - t61 * t80) * MDP(26) + (g(1) * t42 - g(2) * t65 - t62 * t80) * MDP(27);];
taug  = t1;
