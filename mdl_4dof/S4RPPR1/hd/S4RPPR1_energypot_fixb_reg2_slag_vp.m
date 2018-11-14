% Calculate inertial parameters regressor of potential energy for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:47
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4RPPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:46:45
% EndTime: 2018-11-14 13:46:45
% DurationCPUTime: 0.04s
% Computational Cost: add. (62->27), mult. (50->32), div. (0->0), fcn. (44->6), ass. (0->17)
t29 = qJ(2) + pkin(4);
t37 = g(3) * t29;
t28 = qJ(1) + pkin(6);
t24 = sin(t28);
t25 = cos(t28);
t33 = cos(qJ(1));
t36 = t33 * pkin(1) + t25 * pkin(2) + t24 * qJ(3);
t31 = sin(qJ(1));
t35 = t31 * pkin(1) + t24 * pkin(2) - t25 * qJ(3);
t34 = -g(1) * t33 - g(2) * t31;
t32 = cos(qJ(4));
t30 = sin(qJ(4));
t20 = -g(1) * t25 - g(2) * t24;
t19 = g(1) * t24 - g(2) * t25;
t18 = t24 * t32 - t25 * t30;
t17 = -t24 * t30 - t25 * t32;
t1 = [0, 0, 0, 0, 0, 0, t34, g(1) * t31 - g(2) * t33, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, t20, t19, -g(3), t34 * pkin(1) - t37, 0, 0, 0, 0, 0, 0, t20, -g(3), -t19, -g(1) * t36 - g(2) * t35 - t37, 0, 0, 0, 0, 0, 0, g(1) * t17 - g(2) * t18, -g(1) * t18 - g(2) * t17, g(3), -g(1) * (t25 * pkin(3) + t36) - g(2) * (t24 * pkin(3) + t35) - g(3) * (-pkin(5) + t29);];
U_reg  = t1;
