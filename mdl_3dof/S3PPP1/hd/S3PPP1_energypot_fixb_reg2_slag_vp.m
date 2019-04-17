% Calculate inertial parameters regressor of potential energy for
% S3PPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,theta1]';
% 
% Output:
% U_reg [1x(3*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-17 09:48
% Revision: 3acd05283b8979b361f80d69cfa1c98d98241298 (2019-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S3PPP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPP1_energypot_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PPP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPP1_energypot_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-17 09:48:14
% EndTime: 2019-04-17 09:48:14
% DurationCPUTime: 0.03s
% Computational Cost: add. (21->19), mult. (29->16), div. (0->0), fcn. (22->2), ass. (0->8)
t10 = sin(pkin(3));
t11 = cos(pkin(3));
t14 = t11 * pkin(1) + t10 * qJ(2);
t13 = g(3) * qJ(1);
t12 = t10 * pkin(1) - t11 * qJ(2);
t5 = g(1) * t11 + g(2) * t10;
t4 = g(1) * t10 - g(2) * t11;
t1 = [0, 0, 0, 0, 0, 0, -t5, t4, -g(3), -t13, 0, 0, 0, 0, 0, 0, -g(3), t5, -t4, -g(1) * t14 - g(2) * t12 - t13, 0, 0, 0, 0, 0, 0, -g(3), -t4, -t5, -g(1) * (t11 * qJ(3) + t14) - g(2) * (t10 * qJ(3) + t12) - g(3) * (pkin(2) + qJ(1));];
U_reg  = t1;
