% Calculate Gravitation load on the joints for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:51:13
% EndTime: 2019-12-31 21:51:14
% DurationCPUTime: 0.15s
% Computational Cost: add. (256->48), mult. (219->62), div. (0->0), fcn. (174->8), ass. (0->29)
t75 = cos(qJ(3));
t71 = qJ(3) + qJ(4);
t66 = sin(t71);
t68 = cos(t71);
t85 = t68 * pkin(4) + t66 * qJ(5);
t91 = t75 * pkin(3) + t85;
t90 = pkin(2) + t91;
t89 = MDP(19) + MDP(21);
t88 = MDP(20) - MDP(23);
t87 = pkin(4) * t66;
t72 = qJ(1) + qJ(2);
t69 = cos(t72);
t77 = -pkin(8) - pkin(7);
t86 = t69 * t77;
t84 = qJ(5) * t68;
t83 = t90 * t69;
t67 = sin(t72);
t57 = g(1) * t69 + g(2) * t67;
t44 = -g(3) * t68 + t57 * t66;
t82 = t88 * (g(3) * t66 + t57 * t68) + t89 * t44;
t73 = sin(qJ(3));
t81 = -pkin(3) * t73 - t87;
t79 = (-MDP(22) + MDP(6)) * t57 + (t75 * MDP(12) - t73 * MDP(13) - t88 * t66 + t89 * t68 + MDP(5)) * (g(1) * t67 - g(2) * t69);
t78 = (g(1) * t90 + g(2) * t77) * t67;
t76 = cos(qJ(1));
t74 = sin(qJ(1));
t60 = t69 * t84;
t58 = t67 * t84;
t1 = [(g(1) * t74 - g(2) * t76) * MDP(2) + (g(1) * t76 + g(2) * t74) * MDP(3) + (-g(1) * (-t74 * pkin(1) - t86) - g(2) * (t76 * pkin(1) + t83) + t78) * MDP(24) + t79; (g(1) * t86 - g(2) * t83 + t78) * MDP(24) + t79; (-g(3) * t75 + t57 * t73) * MDP(12) + (g(3) * t73 + t57 * t75) * MDP(13) + (-g(1) * (t81 * t69 + t60) - g(2) * (t81 * t67 + t58) - g(3) * t91) * MDP(24) + t82; (-g(1) * (-t69 * t87 + t60) - g(2) * (-t67 * t87 + t58) - g(3) * t85) * MDP(24) + t82; -t44 * MDP(24);];
taug = t1;
