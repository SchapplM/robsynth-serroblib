% Calculate minimal parameter regressor of potential energy for
% S4RPRP3
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
% U_reg [1x15]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:20
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:20:35
% EndTime: 2021-01-15 10:20:35
% DurationCPUTime: 0.04s
% Computational Cost: add. (40->18), mult. (43->25), div. (0->0), fcn. (37->6), ass. (0->15)
t41 = qJ(2) + pkin(4);
t33 = qJ(1) + pkin(6);
t31 = sin(t33);
t32 = cos(t33);
t40 = g(1) * t32 + g(2) * t31;
t36 = sin(qJ(1));
t38 = cos(qJ(1));
t39 = -g(1) * t38 - g(2) * t36;
t37 = cos(qJ(3));
t35 = sin(qJ(3));
t34 = -qJ(4) - pkin(5);
t30 = t37 * pkin(3) + pkin(2);
t29 = -g(3) * t35 - t40 * t37;
t28 = -g(3) * t37 + t40 * t35;
t1 = [0, t39, g(1) * t36 - g(2) * t38, t39 * pkin(1) - g(3) * t41, 0, 0, 0, 0, 0, t29, t28, t29, t28, -g(1) * t31 + g(2) * t32, -g(1) * (t38 * pkin(1) + t32 * t30 - t31 * t34) - g(2) * (t36 * pkin(1) + t31 * t30 + t32 * t34) - g(3) * (t35 * pkin(3) + t41);];
U_reg = t1;
