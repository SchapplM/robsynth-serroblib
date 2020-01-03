% Calculate Gravitation load on the joints for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:35:15
% EndTime: 2019-12-31 21:35:16
% DurationCPUTime: 0.43s
% Computational Cost: add. (173->67), mult. (425->102), div. (0->0), fcn. (457->8), ass. (0->31)
t73 = sin(qJ(3));
t77 = cos(qJ(3));
t79 = cos(qJ(1));
t91 = t79 * t77;
t75 = sin(qJ(1));
t78 = cos(qJ(2));
t93 = t75 * t78;
t65 = t73 * t93 + t91;
t92 = t79 * t73;
t66 = t77 * t93 - t92;
t72 = sin(qJ(5));
t76 = cos(qJ(5));
t107 = t65 * t76 - t66 * t72;
t67 = -t75 * t77 + t78 * t92;
t68 = t75 * t73 + t78 * t91;
t56 = t67 * t76 - t68 * t72;
t57 = t67 * t72 + t68 * t76;
t84 = t72 * t73 + t76 * t77;
t85 = t72 * t77 - t73 * t76;
t86 = t65 * t72 + t66 * t76;
t74 = sin(qJ(2));
t96 = g(3) * t74;
t113 = -(g(1) * t56 + g(2) * t107 - t85 * t96) * MDP(27) + (g(1) * t57 + g(2) * t86 + t84 * t96) * MDP(28);
t109 = MDP(10) - MDP(19);
t105 = MDP(16) + MDP(18);
t104 = MDP(17) - MDP(20);
t88 = g(1) * t79 + g(2) * t75;
t83 = t78 * pkin(2) + t74 * pkin(7) + pkin(1);
t82 = pkin(3) * t77 + qJ(4) * t73 + pkin(2);
t55 = g(1) * t67 + g(2) * t65 + t73 * t96;
t1 = [t88 * MDP(3) + (-g(1) * (-t66 * pkin(3) - t65 * qJ(4)) - g(2) * (t68 * pkin(3) + t67 * qJ(4)) + (-g(1) * pkin(6) - g(2) * t83) * t79 + (-g(2) * pkin(6) + g(1) * t83) * t75) * MDP(21) + (g(1) * t86 - g(2) * t57) * MDP(27) + (g(1) * t107 - g(2) * t56) * MDP(28) - t104 * (g(1) * t65 - g(2) * t67) + t105 * (g(1) * t66 - g(2) * t68) + (t78 * MDP(9) - t109 * t74 + MDP(2)) * (g(1) * t75 - g(2) * t79); ((-t88 * pkin(7) - g(3) * t82) * t78 + (-g(3) * pkin(7) + t88 * t82) * t74) * MDP(21) + t109 * (t88 * t78 + t96) + (-t84 * MDP(27) + t85 * MDP(28) + t104 * t73 - t105 * t77 - MDP(9)) * (g(3) * t78 - t88 * t74); (-g(1) * (-t67 * pkin(3) + t68 * qJ(4)) - g(2) * (-t65 * pkin(3) + t66 * qJ(4)) - (-pkin(3) * t73 + qJ(4) * t77) * t96) * MDP(21) + t104 * (g(1) * t68 + g(2) * t66 + t77 * t96) + t105 * t55 - t113; -t55 * MDP(21); t113;];
taug = t1;
