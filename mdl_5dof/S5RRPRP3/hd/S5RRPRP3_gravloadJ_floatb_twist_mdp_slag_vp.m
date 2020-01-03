% Calculate Gravitation load on the joints for
% S5RRPRP3
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
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:14
% EndTime: 2019-12-31 19:51:14
% DurationCPUTime: 0.13s
% Computational Cost: add. (227->46), mult. (192->55), div. (0->0), fcn. (147->8), ass. (0->26)
t77 = cos(pkin(8));
t74 = pkin(8) + qJ(4);
t69 = sin(t74);
t70 = cos(t74);
t84 = t70 * pkin(4) + t69 * qJ(5);
t96 = -t77 * pkin(3) - pkin(2) - t84;
t95 = MDP(16) + MDP(18);
t94 = MDP(17) - MDP(20);
t75 = qJ(1) + qJ(2);
t71 = sin(t75);
t72 = cos(t75);
t60 = g(1) * t72 + g(2) * t71;
t79 = sin(qJ(1));
t90 = t79 * pkin(1);
t78 = -pkin(7) - qJ(3);
t89 = t72 * t78;
t88 = t72 * pkin(2) + t71 * qJ(3);
t86 = t96 * t72;
t85 = -t71 * pkin(2) + t72 * qJ(3);
t59 = g(1) * t71 - g(2) * t72;
t82 = (-MDP(19) + MDP(6) - MDP(9)) * t60 + (MDP(7) * t77 - MDP(8) * sin(pkin(8)) - t94 * t69 + t95 * t70 + MDP(5)) * t59;
t81 = (-g(1) * t96 + g(2) * t78) * t71;
t80 = cos(qJ(1));
t73 = t80 * pkin(1);
t45 = -g(3) * t70 + t60 * t69;
t1 = [t82 + (-g(1) * (-t89 - t90) - g(2) * (t73 - t86) + t81) * MDP(21) + (-g(1) * (t85 - t90) - g(2) * (t73 + t88)) * MDP(10) + (g(1) * t79 - g(2) * t80) * MDP(2) + (g(1) * t80 + g(2) * t79) * MDP(3); (-g(1) * t85 - g(2) * t88) * MDP(10) + (g(1) * t89 + g(2) * t86 + t81) * MDP(21) + t82; (-MDP(10) - MDP(21)) * t59; (-g(3) * t84 + t60 * (pkin(4) * t69 - qJ(5) * t70)) * MDP(21) + t94 * (g(3) * t69 + t60 * t70) + t95 * t45; -t45 * MDP(21);];
taug = t1;
