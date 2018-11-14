% Calculate minimal parameter regressor of potential energy for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x10]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:54
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:53:29
% EndTime: 2018-11-14 13:53:29
% DurationCPUTime: 0.05s
% Computational Cost: add. (30->14), mult. (19->19), div. (0->0), fcn. (16->6), ass. (0->9)
t68 = qJ(1) + qJ(2);
t70 = cos(qJ(1));
t69 = sin(qJ(1));
t67 = cos(t68);
t66 = sin(t68);
t65 = pkin(7) + qJ(4) + t68;
t64 = cos(t65);
t63 = sin(t65);
t1 = [0, -g(1) * t70 - g(2) * t69, g(1) * t69 - g(2) * t70, 0, -g(1) * t67 - g(2) * t66, g(1) * t66 - g(2) * t67, -g(1) * (t70 * pkin(1) + pkin(2) * t67) - g(2) * (t69 * pkin(1) + pkin(2) * t66) - g(3) * (qJ(3) + pkin(5) + pkin(4)) 0, -g(1) * t64 - g(2) * t63, g(1) * t63 - g(2) * t64;];
U_reg  = t1;
