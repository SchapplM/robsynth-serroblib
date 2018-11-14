% Calculate inertial parameters regressor of potential energy for
% S4PPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:58
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PPRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:57:27
% EndTime: 2018-11-14 13:57:27
% DurationCPUTime: 0.04s
% Computational Cost: add. (40->25), mult. (30->22), div. (0->0), fcn. (20->4), ass. (0->12)
t28 = g(3) * (pkin(4) + qJ(2));
t23 = cos(pkin(5));
t27 = t23 * pkin(2) + pkin(1);
t26 = g(2) * qJ(1);
t22 = sin(pkin(5));
t25 = t22 * pkin(2) + qJ(1);
t21 = pkin(5) + qJ(3);
t18 = cos(t21);
t17 = sin(t21);
t16 = -g(1) * t18 - g(2) * t17;
t15 = g(1) * t17 - g(2) * t18;
t1 = [0, 0, 0, 0, 0, 0, -g(1), g(3), -g(2), -t26, 0, 0, 0, 0, 0, 0, -g(1) * t23 - g(2) * t22, g(1) * t22 - g(2) * t23, -g(3), -g(1) * pkin(1) - g(3) * qJ(2) - t26, 0, 0, 0, 0, 0, 0, t16, t15, -g(3), -g(1) * t27 - g(2) * t25 - t28, 0, 0, 0, 0, 0, 0, t16, -g(3), -t15, -g(1) * (t18 * pkin(3) + t17 * qJ(4) + t27) - g(2) * (t17 * pkin(3) - t18 * qJ(4) + t25) - t28;];
U_reg  = t1;
