% Calculate time derivative of joint inertia matrix for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:17
% EndTime: 2019-12-31 17:32:17
% DurationCPUTime: 0.25s
% Computational Cost: add. (231->59), mult. (677->107), div. (0->0), fcn. (563->6), ass. (0->32)
t21 = sin(pkin(8));
t22 = cos(pkin(8));
t31 = t21 ^ 2 + t22 ^ 2;
t28 = m(5) * t31;
t23 = sin(qJ(5));
t25 = cos(qJ(5));
t12 = t23 * t21 - t25 * t22;
t24 = sin(qJ(3));
t9 = t12 * t24;
t38 = qJD(4) * t28;
t37 = 2 * m(6);
t13 = t25 * t21 + t23 * t22;
t11 = t13 * qJD(5);
t36 = t12 * t11;
t10 = t12 * qJD(5);
t35 = t13 * t10;
t32 = pkin(6) + qJ(4);
t26 = cos(qJ(3));
t29 = qJD(3) * t26;
t27 = t31 * mrSges(5,3);
t5 = t11 * mrSges(6,1) - t10 * mrSges(6,2);
t14 = t32 * t21;
t15 = t32 * t22;
t6 = -t25 * t14 - t23 * t15;
t7 = -t23 * t14 + t25 * t15;
t8 = t13 * t24;
t18 = -t22 * pkin(4) - pkin(3);
t4 = qJD(5) * t9 - t13 * t29;
t3 = -qJD(5) * t8 - t12 * t29;
t2 = -t13 * qJD(4) - t7 * qJD(5);
t1 = -t12 * qJD(4) + t6 * qJD(5);
t16 = [(-t35 + t36) * t37; m(6) * (-t9 * t10 - t8 * t11 + t4 * t12 - t3 * t13); (-t9 * t3 - t8 * t4) * t37 + 0.4e1 * (m(5) * (-0.1e1 + t31) / 0.2e1 - m(6) / 0.2e1) * t24 * t29; m(6) * (-t1 * t13 + t7 * t10 + t6 * t11 + t2 * t12); -t26 * t5 + m(6) * (-t1 * t9 - t2 * t8 + t7 * t3 + t6 * t4) + t24 * t38 + ((-m(5) * pkin(3) + m(6) * t18 - t22 * mrSges(5,1) + t12 * mrSges(6,1) + t21 * mrSges(5,2) + t13 * mrSges(6,2) - mrSges(4,1)) * t24 + (qJ(4) * t28 - mrSges(4,2) + t27) * t26) * qJD(3) + (-t8 * t10 + t9 * t11 - t3 * t12 - t4 * t13) * mrSges(6,3); -0.2e1 * Ifges(6,1) * t35 + 0.2e1 * Ifges(6,2) * t36 + (t7 * t1 + t6 * t2) * t37 + 0.2e1 * t18 * t5 + 0.2e1 * qJ(4) * t38 + 0.2e1 * t27 * qJD(4) + 0.2e1 * (t10 * t12 - t13 * t11) * Ifges(6,4) + 0.2e1 * (-t1 * t12 + t6 * t10 - t7 * t11 - t2 * t13) * mrSges(6,3); 0; (m(5) + m(6)) * t24 * qJD(3); t5; 0; t5; t4 * mrSges(6,1) - t3 * mrSges(6,2); t2 * mrSges(6,1) - t1 * mrSges(6,2) - Ifges(6,5) * t10 - Ifges(6,6) * t11; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t16(1), t16(2), t16(4), t16(7), t16(11); t16(2), t16(3), t16(5), t16(8), t16(12); t16(4), t16(5), t16(6), t16(9), t16(13); t16(7), t16(8), t16(9), t16(10), t16(14); t16(11), t16(12), t16(13), t16(14), t16(15);];
Mq = res;
