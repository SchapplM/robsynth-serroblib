% Calculate minimal parameter regressor of potential energy for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:55
% EndTime: 2024-09-27 21:45:55
% DurationCPUTime: 0.03s
% Computational Cost: add. (139->42), mult. (77->62), div. (0->0), fcn. (76->18), ass. (0->32)
t96 = sin(pkin(5)) * g(3);
t91 = cos(pkin(5));
t92 = sin(qJ(3));
t95 = t91 * t92;
t93 = cos(qJ(3));
t94 = t91 * t93;
t89 = qJ(3) + qJ(4);
t84 = pkin(5) - t89;
t83 = pkin(5) + t89;
t88 = pkin(10) + qJ(2);
t87 = qJ(5) + t89;
t86 = cos(t89);
t85 = sin(t89);
t82 = cos(t88);
t81 = sin(t88);
t80 = cos(t87);
t79 = sin(t87);
t78 = -qJ(5) + t84;
t77 = qJ(5) + t83;
t76 = cos(t83);
t75 = sin(t84);
t74 = cos(t77);
t73 = sin(t78);
t72 = cos(t84) / 0.2e1;
t71 = sin(t83) / 0.2e1;
t70 = cos(t78) / 0.2e1;
t69 = sin(t77) / 0.2e1;
t68 = t72 + t76 / 0.2e1;
t67 = t71 - t75 / 0.2e1;
t66 = t70 + t74 / 0.2e1;
t65 = t69 - t73 / 0.2e1;
t1 = [-g(3) * qJ(1), 0, -g(1) * t82 - g(2) * t81, g(1) * t81 - g(2) * t82, 0, 0, 0, 0, 0, -g(1) * (-t81 * t95 + t82 * t93) - g(2) * (t81 * t93 + t82 * t95) - t92 * t96, -g(1) * (-t81 * t94 - t82 * t92) - g(2) * (-t81 * t92 + t82 * t94) - t93 * t96, 0, 0, 0, 0, 0, -g(1) * (-t81 * t67 + t82 * t86) - g(2) * (t82 * t67 + t81 * t86) - g(3) * (t72 - t76 / 0.2e1), -g(1) * (-t81 * t68 - t82 * t85) - g(2) * (t82 * t68 - t81 * t85) - g(3) * (t71 + t75 / 0.2e1), 0, 0, 0, 0, 0, -g(1) * (-t81 * t65 + t82 * t80) - g(2) * (t82 * t65 + t81 * t80) - g(3) * (t70 - t74 / 0.2e1), -g(1) * (-t81 * t66 - t82 * t79) - g(2) * (t82 * t66 - t81 * t79) - g(3) * (t69 + t73 / 0.2e1);];
U_reg = t1;
