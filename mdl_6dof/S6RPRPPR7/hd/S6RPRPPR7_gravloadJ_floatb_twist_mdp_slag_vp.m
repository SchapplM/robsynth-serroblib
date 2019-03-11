% Calculate Gravitation load on the joints for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:57:36
% EndTime: 2019-03-09 02:57:37
% DurationCPUTime: 0.22s
% Computational Cost: add. (142->59), mult. (210->78), div. (0->0), fcn. (175->8), ass. (0->30)
t62 = sin(qJ(1));
t65 = cos(qJ(1));
t83 = -g(1) * t62 + g(2) * t65;
t58 = qJ(3) + pkin(9);
t52 = sin(t58);
t80 = g(3) * t52;
t61 = sin(qJ(3));
t79 = t61 * pkin(3);
t60 = sin(qJ(6));
t78 = t62 * t60;
t63 = cos(qJ(6));
t77 = t62 * t63;
t76 = t65 * t60;
t75 = t65 * t63;
t74 = t65 * pkin(1) + t62 * qJ(2);
t73 = -MDP(15) - MDP(19);
t72 = -t62 * pkin(1) + t65 * qJ(2);
t48 = g(1) * t65 + g(2) * t62;
t53 = cos(t58);
t71 = t52 * pkin(4) - t53 * qJ(5);
t59 = -qJ(4) - pkin(7);
t70 = t62 * t59 + t65 * t79 + t72;
t69 = -t65 * t59 + t62 * t79 + t74;
t64 = cos(qJ(3));
t46 = -t53 * t78 + t75;
t45 = -t53 * t77 - t76;
t44 = -t53 * t76 - t77;
t43 = -t53 * t75 + t78;
t42 = -t53 * t83 - t80;
t1 = [(-g(1) * t72 - g(2) * t74) * MDP(6) + (-g(1) * t70 - g(2) * t69) * MDP(15) + (-g(1) * (t71 * t65 + t70) - g(2) * (t71 * t62 + t69)) * MDP(19) + (-g(1) * t44 - g(2) * t46) * MDP(25) + (-g(1) * t43 - g(2) * t45) * MDP(26) - (MDP(2) - MDP(4) + MDP(14) + MDP(16)) * t83 + (-t61 * MDP(12) - t64 * MDP(13) + t52 * MDP(17) + t53 * MDP(18) + MDP(3) - MDP(5)) * t48; -(-MDP(6) + t73) * t83; (g(3) * t64 - t61 * t83) * MDP(13) + t42 * MDP(17) + (-g(3) * (-t71 - t79) + t83 * (pkin(3) * t64 + pkin(4) * t53 + qJ(5) * t52)) * MDP(19) + (pkin(3) * MDP(15) + MDP(12)) * (g(3) * t61 + t83 * t64) + (MDP(25) * t60 + MDP(26) * t63 + MDP(18)) * (-g(3) * t53 + t52 * t83); t73 * t48; t42 * MDP(19); (-g(1) * t45 + g(2) * t43 - t63 * t80) * MDP(25) + (g(1) * t46 - g(2) * t44 + t60 * t80) * MDP(26);];
taug  = t1;
