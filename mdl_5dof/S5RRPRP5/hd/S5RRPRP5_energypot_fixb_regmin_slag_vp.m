% Calculate minimal parameter regressor of potential energy for
% S5RRPRP5
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
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:19
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:18:17
% EndTime: 2021-01-15 20:18:17
% DurationCPUTime: 0.06s
% Computational Cost: add. (94->32), mult. (81->38), div. (0->0), fcn. (75->8), ass. (0->20)
t70 = sin(qJ(2));
t76 = t70 * pkin(2) + pkin(5);
t72 = cos(qJ(2));
t60 = t72 * pkin(2) + pkin(1);
t69 = pkin(6) + qJ(3);
t68 = qJ(2) + pkin(8);
t71 = sin(qJ(1));
t73 = cos(qJ(1));
t75 = g(1) * t73 + g(2) * t71;
t63 = qJ(4) + t68;
t58 = sin(t63);
t59 = cos(t63);
t62 = cos(t68);
t74 = pkin(3) * t62 + pkin(4) * t59 + qJ(5) * t58 + t60;
t67 = -pkin(7) - t69;
t61 = sin(t68);
t57 = g(1) * t71 - g(2) * t73;
t55 = -g(3) * t58 - t75 * t59;
t54 = -g(3) * t59 + t75 * t58;
t1 = [0, -t75, t57, 0, 0, 0, 0, 0, -g(3) * t70 - t75 * t72, -g(3) * t72 + t75 * t70, -g(3) * t61 - t75 * t62, -g(3) * t62 + t75 * t61, -t57, -g(1) * (t73 * t60 + t69 * t71) - g(2) * (t71 * t60 - t73 * t69) - g(3) * t76, 0, 0, 0, 0, 0, t55, t54, t55, -t57, -t54, -g(3) * (pkin(3) * t61 + t58 * pkin(4) - t59 * qJ(5) + t76) + (-g(1) * t74 - g(2) * t67) * t73 + (g(1) * t67 - g(2) * t74) * t71;];
U_reg = t1;
