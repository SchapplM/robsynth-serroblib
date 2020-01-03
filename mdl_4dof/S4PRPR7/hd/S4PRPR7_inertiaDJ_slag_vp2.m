% Calculate time derivative of joint inertia matrix for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR7_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR7_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:37
% EndTime: 2019-12-31 16:25:38
% DurationCPUTime: 0.11s
% Computational Cost: add. (45->29), mult. (141->48), div. (0->0), fcn. (66->4), ass. (0->16)
t19 = m(4) + m(5);
t7 = sin(qJ(4));
t9 = cos(qJ(4));
t2 = t7 * mrSges(5,1) + t9 * mrSges(5,2);
t18 = mrSges(4,3) + t2;
t6 = t9 ^ 2;
t16 = t7 ^ 2 + t6;
t8 = sin(qJ(2));
t15 = qJD(2) * t8;
t10 = cos(qJ(2));
t14 = qJD(2) * t10;
t13 = qJD(4) * t10;
t12 = m(5) * t16;
t11 = -pkin(2) - pkin(5);
t1 = (mrSges(5,1) * t9 - mrSges(5,2) * t7) * qJD(4);
t3 = [0.2e1 * m(5) * (0.1e1 - t16) * t8 * t14; t8 * t1 + ((-mrSges(3,2) + t18) * t10 + (-m(4) * pkin(2) - t16 * mrSges(5,3) + t11 * t12 - mrSges(3,1) + mrSges(4,2)) * t8) * qJD(2) + t19 * (qJ(3) * t14 + t8 * qJD(3)); 0.2e1 * qJ(3) * t1 + 0.2e1 * (t19 * qJ(3) + t18) * qJD(3) + 0.2e1 * (-t6 * Ifges(5,4) + (Ifges(5,4) * t7 + (-Ifges(5,1) + Ifges(5,2)) * t9) * t7) * qJD(4); (m(4) + t12) * t15; 0; 0; (t9 * t13 - t7 * t15) * mrSges(5,2) + (t7 * t13 + t9 * t15) * mrSges(5,1); (-Ifges(5,5) * t7 - Ifges(5,6) * t9 - t2 * t11) * qJD(4); -t2 * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1), t3(2), t3(4), t3(7); t3(2), t3(3), t3(5), t3(8); t3(4), t3(5), t3(6), t3(9); t3(7), t3(8), t3(9), t3(10);];
Mq = res;
