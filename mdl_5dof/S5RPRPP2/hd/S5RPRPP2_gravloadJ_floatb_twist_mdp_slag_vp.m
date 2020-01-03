% Calculate Gravitation load on the joints for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPRPP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:12
% EndTime: 2019-12-31 18:11:13
% DurationCPUTime: 0.22s
% Computational Cost: add. (149->46), mult. (174->61), div. (0->0), fcn. (135->6), ass. (0->27)
t79 = MDP(10) + MDP(12) + MDP(16);
t78 = MDP(11) - MDP(14) - MDP(17);
t57 = qJ(1) + pkin(7);
t51 = sin(t57);
t52 = cos(t57);
t43 = g(1) * t52 + g(2) * t51;
t58 = sin(qJ(3));
t77 = t43 * t58;
t53 = t58 * qJ(4);
t60 = cos(qJ(3));
t69 = t60 * pkin(3) + t53;
t75 = pkin(3) * t58;
t74 = g(1) * t51;
t71 = t60 * pkin(4);
t70 = t52 * t60;
t68 = qJ(4) * t60;
t67 = -MDP(15) - MDP(19);
t59 = sin(qJ(1));
t66 = -t59 * pkin(1) + t52 * pkin(6);
t61 = cos(qJ(1));
t65 = t61 * pkin(1) + pkin(3) * t70 + t51 * pkin(6) + (pkin(2) + t53) * t52;
t64 = -g(2) * t52 + t74;
t62 = -pkin(2) - t69;
t46 = t52 * t68;
t44 = t51 * t68;
t39 = -g(3) * t60 + t77;
t1 = [(g(1) * t61 + g(2) * t59) * MDP(3) + (-g(1) * t66 - g(2) * t65 - t62 * t74) * MDP(15) + (-g(1) * (-t52 * qJ(5) + t66) - g(2) * (pkin(4) * t70 + t65) + (-g(1) * (t62 - t71) + g(2) * qJ(5)) * t51) * MDP(19) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t59 - g(2) * t61) + (-MDP(13) + MDP(18)) * t43 + (-t78 * t58 + t79 * t60) * t64; (-MDP(4) + t67) * g(3); (-g(1) * (-t52 * t75 + t46) - g(2) * (-t51 * t75 + t44) - g(3) * t69) * MDP(15) + (-g(1) * t46 - g(2) * t44 - g(3) * (t69 + t71) + (pkin(3) + pkin(4)) * t77) * MDP(19) + t78 * (g(3) * t58 + t43 * t60) + t79 * t39; t67 * t39; t64 * MDP(19);];
taug = t1;
