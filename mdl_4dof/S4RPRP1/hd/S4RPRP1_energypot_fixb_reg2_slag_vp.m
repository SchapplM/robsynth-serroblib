% Calculate inertial parameters regressor of potential energy for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4RPRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:48:38
% EndTime: 2018-11-14 13:48:38
% DurationCPUTime: 0.04s
% Computational Cost: add. (64->27), mult. (40->28), div. (0->0), fcn. (30->6), ass. (0->16)
t38 = qJ(2) + pkin(4);
t39 = g(3) * (pkin(5) + t38);
t32 = qJ(1) + pkin(6);
t26 = sin(t32);
t33 = sin(qJ(1));
t37 = t33 * pkin(1) + pkin(2) * t26;
t27 = cos(t32);
t34 = cos(qJ(1));
t36 = t34 * pkin(1) + pkin(2) * t27;
t35 = -g(1) * t34 - g(2) * t33;
t28 = qJ(3) + t32;
t25 = cos(t28);
t24 = sin(t28);
t21 = -g(1) * t25 - g(2) * t24;
t20 = g(1) * t24 - g(2) * t25;
t1 = [0, 0, 0, 0, 0, 0, t35, g(1) * t33 - g(2) * t34, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, -g(1) * t27 - g(2) * t26, g(1) * t26 - g(2) * t27, -g(3), t35 * pkin(1) - g(3) * t38, 0, 0, 0, 0, 0, 0, t21, t20, -g(3), -g(1) * t36 - g(2) * t37 - t39, 0, 0, 0, 0, 0, 0, t21, -g(3), -t20, -g(1) * (t25 * pkin(3) + t24 * qJ(4) + t36) - g(2) * (t24 * pkin(3) - t25 * qJ(4) + t37) - t39;];
U_reg  = t1;
