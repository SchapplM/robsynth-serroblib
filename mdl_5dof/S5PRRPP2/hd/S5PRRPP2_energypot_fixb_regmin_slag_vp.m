% Calculate minimal parameter regressor of potential energy for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:33
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:32:31
% EndTime: 2021-01-15 15:32:31
% DurationCPUTime: 0.10s
% Computational Cost: add. (109->45), mult. (149->64), div. (0->0), fcn. (158->8), ass. (0->30)
t66 = qJ(4) + pkin(6);
t68 = sin(qJ(2));
t83 = t66 * t68 + pkin(1);
t82 = g(3) * t68;
t64 = sin(pkin(7));
t67 = sin(qJ(3));
t81 = t64 * t67;
t70 = cos(qJ(2));
t80 = t64 * t70;
t65 = cos(pkin(7));
t79 = t65 * t70;
t77 = t67 * t70;
t69 = cos(qJ(3));
t76 = t69 * t70;
t56 = t69 * pkin(3) + pkin(2);
t75 = t68 * t56 - t70 * t66 + qJ(1);
t74 = pkin(3) * t81 + t64 * pkin(5) + t56 * t79 + t83 * t65;
t73 = g(1) * t65 + g(2) * t64;
t63 = qJ(3) + pkin(8);
t57 = sin(t63);
t58 = cos(t63);
t43 = t57 * t80 + t65 * t58;
t45 = t57 * t79 - t64 * t58;
t72 = g(1) * t45 + g(2) * t43 + t57 * t82;
t71 = t56 * t80 + (-pkin(3) * t67 - pkin(5)) * t65 + t83 * t64;
t47 = -g(3) * t70 + t73 * t68;
t46 = t64 * t57 + t58 * t79;
t44 = -t65 * t57 + t58 * t80;
t42 = -g(1) * t46 - g(2) * t44 - t58 * t82;
t1 = [-g(3) * qJ(1), 0, -t73 * t70 - t82, t47, 0, 0, 0, 0, 0, -g(1) * (t65 * t76 + t81) - g(2) * (t64 * t76 - t65 * t67) - t69 * t82, -g(1) * (t64 * t69 - t65 * t77) - g(2) * (-t64 * t77 - t65 * t69) + t67 * t82, t42, t72, -t47, -g(1) * t74 - g(2) * t71 - g(3) * t75, t42, -t47, -t72, -g(1) * (t46 * pkin(4) + t45 * qJ(5) + t74) - g(2) * (t44 * pkin(4) + t43 * qJ(5) + t71) - g(3) * ((pkin(4) * t58 + qJ(5) * t57) * t68 + t75);];
U_reg = t1;
