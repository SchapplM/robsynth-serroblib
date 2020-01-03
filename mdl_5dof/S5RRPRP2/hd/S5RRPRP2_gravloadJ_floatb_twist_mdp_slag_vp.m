% Calculate Gravitation load on the joints for
% S5RRPRP2
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
%   see S5RRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:43
% EndTime: 2019-12-31 19:49:44
% DurationCPUTime: 0.14s
% Computational Cost: add. (194->42), mult. (158->54), div. (0->0), fcn. (119->8), ass. (0->27)
t65 = sin(qJ(4));
t67 = cos(qJ(4));
t72 = t67 * pkin(4) + t65 * qJ(5);
t86 = -pkin(3) - t72;
t85 = MDP(13) + MDP(15);
t84 = MDP(14) - MDP(17);
t64 = qJ(1) + qJ(2);
t60 = pkin(8) + t64;
t57 = sin(t60);
t58 = cos(t60);
t75 = g(1) * t58 + g(2) * t57;
t61 = sin(t64);
t83 = pkin(2) * t61;
t82 = g(1) * t57;
t62 = cos(t64);
t59 = pkin(2) * t62;
t77 = t57 * pkin(7) - t86 * t58 + t59;
t66 = sin(qJ(1));
t76 = -t66 * pkin(1) - t83;
t73 = g(1) * t61 - g(2) * t62;
t70 = -t75 * MDP(16) + t73 * MDP(5) + (g(1) * t62 + g(2) * t61) * MDP(6) + (-t84 * t65 + t85 * t67) * (-g(2) * t58 + t82);
t69 = t86 * t82;
t68 = cos(qJ(1));
t63 = t68 * pkin(1);
t55 = t58 * pkin(7);
t41 = -g(3) * t67 + t65 * t75;
t1 = [(g(1) * t66 - g(2) * t68) * MDP(2) + (g(1) * t68 + g(2) * t66) * MDP(3) + (-g(1) * t76 - g(2) * (t59 + t63)) * MDP(7) + (-g(1) * (t55 + t76) - g(2) * (t63 + t77) - t69) * MDP(18) + t70; t73 * pkin(2) * MDP(7) + (-g(1) * (t55 - t83) - g(2) * t77 - t69) * MDP(18) + t70; (-MDP(18) - MDP(7)) * g(3); (-g(3) * t72 + t75 * (pkin(4) * t65 - qJ(5) * t67)) * MDP(18) + t84 * (g(3) * t65 + t67 * t75) + t85 * t41; -t41 * MDP(18);];
taug = t1;
