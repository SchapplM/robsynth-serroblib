% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:00
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PPRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:59:33
% EndTime: 2018-11-14 13:59:33
% DurationCPUTime: 0.07s
% Computational Cost: add. (31->14), mult. (92->35), div. (0->0), fcn. (58->6), ass. (0->15)
t10 = cos(pkin(6));
t12 = sin(qJ(3));
t14 = cos(qJ(3));
t9 = sin(pkin(6));
t4 = (t10 * t14 - t12 * t9) * qJD(1);
t15 = qJD(1) ^ 2;
t13 = cos(qJ(4));
t11 = sin(qJ(4));
t8 = qJD(2) ^ 2 / 0.2e1;
t7 = qJD(3) + qJD(4);
t5 = (t10 * t12 + t14 * t9) * qJD(1);
t3 = qJD(3) * pkin(3) + t4;
t2 = t11 * t3 + t13 * t5;
t1 = -t11 * t5 + t13 * t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t15 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 + (t9 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1) * t15, 0, 0, 0, 0, 0, qJD(3) ^ 2 / 0.2e1, t4 * qJD(3), -t5 * qJD(3), 0, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t8, 0, 0, 0, 0, 0, t7 ^ 2 / 0.2e1, t1 * t7, -t2 * t7, 0, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t8;];
T_reg  = t6;
