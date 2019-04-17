% Calculate minimal parameter regressor of potential energy for
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
% U_reg [1x3]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-17 09:48
% Revision: 3acd05283b8979b361f80d69cfa1c98d98241298 (2019-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S3PPP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPP1_energypot_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PPP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPP1_energypot_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-17 09:48:14
% EndTime: 2019-04-17 09:48:14
% DurationCPUTime: 0.02s
% Computational Cost: add. (12->10), mult. (17->12), div. (0->0), fcn. (10->2), ass. (0->6)
t25 = g(3) * qJ(1);
t21 = sin(pkin(3));
t22 = cos(pkin(3));
t24 = t22 * pkin(1) + t21 * qJ(2);
t23 = t21 * pkin(1) - t22 * qJ(2);
t1 = [-t25, -g(1) * t24 - g(2) * t23 - t25, -g(1) * (t22 * qJ(3) + t24) - g(2) * (t21 * qJ(3) + t23) - g(3) * (pkin(2) + qJ(1));];
U_reg  = t1;
