% Calculate inertial parameters regressor of potential energy for
% S3PPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% 
% Output:
% U_reg [1x(3*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:10
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S3PPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR1_energypot_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR1_energypot_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:10:05
% EndTime: 2018-11-14 10:10:05
% DurationCPUTime: 0.02s
% Computational Cost: add. (13->13), mult. (11->9), div. (0->0), fcn. (4->2), ass. (0->5)
t9 = g(1) * qJ(2);
t8 = g(2) * qJ(1);
t7 = cos(qJ(3));
t6 = sin(qJ(3));
t1 = [0, 0, 0, 0, 0, 0, g(3), g(1), -g(2), -t8, 0, 0, 0, 0, 0, 0, -g(2), -g(3), -g(1), g(3) * pkin(1) - t8 - t9, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t7, -g(1) * t7 + g(2) * t6, g(3), -t9 - g(2) * (pkin(2) + qJ(1)) - g(3) * (-pkin(3) - pkin(1));];
U_reg  = t1;
