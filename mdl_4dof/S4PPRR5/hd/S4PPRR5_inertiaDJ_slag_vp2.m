% Calculate time derivative of joint inertia matrix for
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR5_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR5_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:41
% EndTime: 2019-12-31 16:19:42
% DurationCPUTime: 0.13s
% Computational Cost: add. (47->25), mult. (188->59), div. (0->0), fcn. (94->4), ass. (0->19)
t10 = cos(qJ(4));
t8 = sin(qJ(4));
t4 = t8 ^ 2;
t18 = t10 ^ 2 + t4;
t26 = pkin(5) * t18;
t2 = (mrSges(5,1) * t8 + mrSges(5,2) * t10) * qJD(4);
t25 = -t2 + (t18 * mrSges(5,3) - mrSges(4,2)) * qJD(3);
t9 = sin(qJ(3));
t24 = m(5) * t9;
t23 = -0.1e1 + t18;
t3 = -t10 * mrSges(5,1) + t8 * mrSges(5,2);
t21 = -m(5) * pkin(3) - mrSges(4,1) + t3;
t11 = cos(qJ(3));
t14 = qJD(3) * t11;
t20 = t23 * t14 * t24;
t15 = qJD(3) * t9;
t13 = qJD(4) * t10;
t12 = qJD(4) * t11;
t1 = [-0.2e1 * t20; m(5) * t23 * (t11 ^ 2 - t9 ^ 2) * qJD(3); 0.2e1 * t20; (t21 * t11 - t24 * t26) * qJD(3) - t25 * t9; (m(5) * t11 * t26 + t21 * t9) * qJD(3) + t25 * t11; -0.2e1 * t4 * Ifges(5,4) * qJD(4) - 0.2e1 * pkin(3) * t2 + 0.2e1 * (Ifges(5,4) * t10 + (Ifges(5,1) - Ifges(5,2)) * t8) * t13; (t10 * t15 + t8 * t12) * mrSges(5,2) + (-t10 * t12 + t8 * t15) * mrSges(5,1); (qJD(4) * t8 * t9 - t10 * t14) * mrSges(5,2) + (-t9 * t13 - t8 * t14) * mrSges(5,1); (Ifges(5,5) * t10 - Ifges(5,6) * t8 + t3 * pkin(5)) * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
