% Calculate inertial parameters regressor of coriolis matrix for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPP1_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:51:42
% EndTime: 2018-11-14 13:51:42
% DurationCPUTime: 0.21s
% Computational Cost: add. (134->45), mult. (315->56), div. (0->0), fcn. (237->4), ass. (0->36)
t22 = sin(pkin(6));
t42 = t22 * pkin(2);
t25 = cos(qJ(2));
t41 = t25 * pkin(1);
t20 = pkin(2) + t41;
t40 = t22 * t20;
t39 = t22 * t25;
t23 = cos(pkin(6));
t24 = sin(qJ(2));
t38 = t23 * t24;
t37 = pkin(1) * qJD(1);
t36 = pkin(1) * qJD(2);
t13 = (t38 + t39) * pkin(1);
t17 = t22 * t24 * pkin(1);
t14 = t23 * t41 - t17;
t29 = t23 * t20 - t17;
t26 = pkin(1) * t38 + t40;
t6 = qJ(4) + t26;
t1 = t6 * t14 + (-pkin(3) - t29) * t13;
t35 = t1 * qJD(1);
t2 = -t29 * t13 + t26 * t14;
t34 = t2 * qJD(1);
t33 = t6 * qJD(1);
t30 = t14 * qJD(2);
t32 = t30 + qJD(4);
t7 = t13 * qJD(1);
t31 = t14 * qJD(1);
t21 = qJD(1) + qJD(2);
t28 = pkin(1) * t21;
t18 = qJ(4) + t42;
t5 = -qJ(4) + (t41 / 0.2e1 - pkin(2) / 0.2e1 - t20 / 0.2e1) * t22;
t27 = t5 * qJD(1) - t18 * qJD(2);
t8 = t13 * qJD(2);
t4 = t42 / 0.2e1 + qJ(4) + t40 / 0.2e1 + (t38 + t39 / 0.2e1) * pkin(1);
t3 = -t8 - t7;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24 * t36, -t25 * t36, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t30, 0, t2 * qJD(2), 0, 0, 0, 0, 0, 0, -t8, 0, t32, t1 * qJD(2) + t6 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24 * t28, -t25 * t28, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t31 - t30, 0, t34 + (-t13 * t23 + t14 * t22) * qJD(2) * pkin(2), 0, 0, 0, 0, 0, 0, t3, 0, t31 + t32, t35 + (t14 * t18 + t13 * (-t23 * pkin(2) - pkin(3))) * qJD(2) + t4 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t4 * qJD(2) + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t37, t25 * t37, 0, 0, 0, 0, 0, 0, 0, 0, t7, t31, 0, -t34, 0, 0, 0, 0, 0, 0, t7, 0, -t31 + qJD(4), -t5 * qJD(4) - t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t18 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t5 * qJD(2) - t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t9;
