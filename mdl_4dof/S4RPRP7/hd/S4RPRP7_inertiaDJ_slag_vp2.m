% Calculate time derivative of joint inertia matrix for
% S4RPRP7
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_inertiaDJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP7_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP7_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:03
% EndTime: 2019-12-31 16:47:03
% DurationCPUTime: 0.13s
% Computational Cost: add. (66->36), mult. (152->49), div. (0->0), fcn. (66->2), ass. (0->14)
t3 = sin(qJ(3));
t4 = cos(qJ(3));
t16 = t3 * mrSges(4,1) + t4 * mrSges(4,2);
t7 = -t3 * pkin(3) + t4 * qJ(4);
t2 = qJ(2) - t7;
t15 = 0.2e1 * t2;
t14 = m(5) * t3;
t11 = Ifges(5,5) - Ifges(4,4);
t10 = 0.2e1 * t4;
t5 = -pkin(1) - pkin(5);
t9 = (m(5) * t5 - mrSges(5,2)) * t3;
t8 = -t3 * mrSges(5,1) + t4 * mrSges(5,3);
t6 = m(5) * t7 - t16 + t8;
t1 = [0.2e1 * (m(5) * t2 - t8) * (-t4 * qJD(4) + qJD(2) + (pkin(3) * t4 + qJ(4) * t3) * qJD(3)) + 0.2e1 * (mrSges(3,3) + (m(3) + m(4)) * qJ(2) + t16) * qJD(2) + ((0.2e1 * qJ(2) * mrSges(4,1) + mrSges(5,1) * t15 + t11 * t10) * t4 + (-0.2e1 * qJ(2) * mrSges(4,2) + mrSges(5,3) * t15 - 0.2e1 * t11 * t3 + (-Ifges(4,1) - Ifges(5,1) + Ifges(4,2) + Ifges(5,3)) * t10) * t3) * qJD(3); 0; 0; qJD(4) * t9 + ((-qJ(4) * mrSges(5,2) - Ifges(4,6) + Ifges(5,6)) * t4 + (pkin(3) * mrSges(5,2) - Ifges(5,4) - Ifges(4,5)) * t3 + t6 * t5) * qJD(3); t6 * qJD(3) + qJD(4) * t14; 0.2e1 * (m(5) * qJ(4) + mrSges(5,3)) * qJD(4); qJD(3) * t9; qJD(3) * t14; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
