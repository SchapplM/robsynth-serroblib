% Calculate time derivative of joint inertia matrix for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_inertiaDJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP6_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP6_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:59
% EndTime: 2019-12-31 16:46:00
% DurationCPUTime: 0.14s
% Computational Cost: add. (97->46), mult. (203->67), div. (0->0), fcn. (98->2), ass. (0->20)
t9 = cos(qJ(3));
t13 = 0.2e1 * t9;
t19 = 2 * mrSges(5,3);
t18 = 2 * mrSges(4,1);
t17 = Ifges(5,4) + Ifges(4,4);
t8 = sin(qJ(3));
t16 = qJD(3) * t8;
t15 = qJD(3) * t9;
t10 = -pkin(1) - pkin(5);
t14 = qJ(4) - t10;
t12 = m(5) * pkin(3) + mrSges(5,1);
t4 = t14 * t9;
t1 = -t9 * qJD(4) + t14 * t16;
t2 = -qJD(3) * t4 - t8 * qJD(4);
t11 = -t9 * t1 - t8 * t2;
t7 = mrSges(5,1) * t15;
t6 = t8 * pkin(3) + qJ(2);
t5 = pkin(3) * t15 + qJD(2);
t3 = t14 * t8;
t20 = [0.2e1 * t5 * (t8 * mrSges(5,1) + t9 * mrSges(5,2)) + 0.2e1 * t6 * t7 + 0.2e1 * m(5) * (-t4 * t1 - t3 * t2 + t6 * t5) + t11 * t19 + (t8 * t18 + mrSges(4,2) * t13 + (2 * mrSges(3,3)) + 0.2e1 * (m(3) + m(4)) * qJ(2)) * qJD(2) + ((qJ(2) * t18 - t17 * t13 + t3 * t19) * t9 + (-0.2e1 * qJ(2) * mrSges(4,2) - 0.2e1 * t6 * mrSges(5,2) - 0.2e1 * t4 * mrSges(5,3) + 0.2e1 * t17 * t8 + (-Ifges(4,1) - Ifges(5,1) + Ifges(4,2) + Ifges(5,2)) * t13) * t8) * qJD(3); m(5) * ((-t3 * t9 + t4 * t8) * qJD(3) - t11); 0; -t2 * mrSges(5,2) + t12 * t1 + ((-mrSges(4,2) * t10 - Ifges(4,6) - Ifges(5,6)) * t9 + (-mrSges(4,1) * t10 + mrSges(5,3) * pkin(3) - Ifges(4,5) - Ifges(5,5)) * t8) * qJD(3); ((-mrSges(4,2) - mrSges(5,2)) * t9 + (-mrSges(4,1) - t12) * t8) * qJD(3); 0; m(5) * t5 - mrSges(5,2) * t16 + t7; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t20(1), t20(2), t20(4), t20(7); t20(2), t20(3), t20(5), t20(8); t20(4), t20(5), t20(6), t20(9); t20(7), t20(8), t20(9), t20(10);];
Mq = res;
