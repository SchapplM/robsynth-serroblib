% Calculate minimal parameter regressor of potential energy for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% 
% Output:
% U_reg [1x9]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S3RPP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_energypot_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_energypot_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:13:43
% EndTime: 2018-11-14 10:13:43
% DurationCPUTime: 0.02s
% Computational Cost: add. (17->15), mult. (28->16), div. (0->0), fcn. (22->2), ass. (0->7)
t28 = sin(qJ(1));
t29 = cos(qJ(1));
t31 = t29 * pkin(1) + t28 * qJ(2);
t30 = t28 * pkin(1) - t29 * qJ(2);
t23 = g(1) * t29 + g(2) * t28;
t22 = g(1) * t28 - g(2) * t29;
t1 = [0, -t23, t22, t23, -t22, -g(3) * pkin(3) - g(1) * t31 - g(2) * t30, -t22, -t23, -g(1) * (t29 * qJ(3) + t31) - g(2) * (t28 * qJ(3) + t30) - g(3) * (pkin(2) + pkin(3));];
U_reg  = t1;
