% Calculate time derivative of joint inertia matrix for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:25
% EndTime: 2019-12-31 17:33:25
% DurationCPUTime: 0.11s
% Computational Cost: add. (46->29), mult. (145->48), div. (0->0), fcn. (68->4), ass. (0->16)
t19 = m(5) + m(6);
t7 = sin(qJ(5));
t9 = cos(qJ(5));
t2 = t7 * mrSges(6,1) + t9 * mrSges(6,2);
t18 = mrSges(5,3) + t2;
t6 = t9 ^ 2;
t16 = t7 ^ 2 + t6;
t8 = sin(qJ(3));
t15 = qJD(3) * t8;
t10 = cos(qJ(3));
t14 = qJD(3) * t10;
t13 = qJD(5) * t10;
t12 = m(6) * t16;
t11 = -pkin(3) - pkin(6);
t1 = (mrSges(6,1) * t9 - mrSges(6,2) * t7) * qJD(5);
t3 = [0; 0; 0.2e1 * m(6) * (0.1e1 - t16) * t8 * t14; 0; t8 * t1 + ((-mrSges(4,2) + t18) * t10 + (-m(5) * pkin(3) - t16 * mrSges(6,3) + t11 * t12 - mrSges(4,1) + mrSges(5,2)) * t8) * qJD(3) + t19 * (qJ(4) * t14 + qJD(4) * t8); 0.2e1 * qJ(4) * t1 + 0.2e1 * (t19 * qJ(4) + t18) * qJD(4) + 0.2e1 * (-t6 * Ifges(6,4) + (Ifges(6,4) * t7 + (-Ifges(6,1) + Ifges(6,2)) * t9) * t7) * qJD(5); 0; (m(5) + t12) * t15; 0; 0; t1; (t9 * t13 - t7 * t15) * mrSges(6,2) + (t7 * t13 + t9 * t15) * mrSges(6,1); (-Ifges(6,5) * t7 - Ifges(6,6) * t9 - t2 * t11) * qJD(5); -t2 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
