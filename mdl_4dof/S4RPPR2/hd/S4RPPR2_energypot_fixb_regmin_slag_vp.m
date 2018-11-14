% Calculate minimal parameter regressor of potential energy for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% U_reg [1x12]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4RPPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:47:31
% EndTime: 2018-11-14 13:47:32
% DurationCPUTime: 0.06s
% Computational Cost: add. (35->21), mult. (48->32), div. (0->0), fcn. (50->6), ass. (0->16)
t68 = sin(qJ(1));
t69 = cos(qJ(1));
t71 = pkin(1) * t69 + qJ(2) * t68;
t70 = pkin(1) * t68 - qJ(2) * t69;
t67 = cos(pkin(6));
t66 = sin(pkin(6));
t65 = pkin(6) + qJ(4);
t61 = cos(t65);
t60 = sin(t65);
t59 = -g(1) * t69 - g(2) * t68;
t58 = g(1) * t68 - g(2) * t69;
t57 = -t66 * t69 + t67 * t68;
t56 = -t66 * t68 - t67 * t69;
t55 = -t60 * t69 + t61 * t68;
t54 = -t60 * t68 - t61 * t69;
t1 = [0, t59, t58, t59, -t58, -g(3) * pkin(4) - g(1) * t71 - g(2) * t70, g(1) * t56 - g(2) * t57, -g(1) * t57 - g(2) * t56, -g(1) * (pkin(2) * t69 + t71) - g(2) * (pkin(2) * t68 + t70) - g(3) * (-qJ(3) + pkin(4)) 0, g(1) * t54 - g(2) * t55, -g(1) * t55 - g(2) * t54;];
U_reg  = t1;
