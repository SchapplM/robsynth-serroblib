% Calculate time derivative of joint inertia matrix for
% S4PRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
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
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_inertiaDJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP2_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP2_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP2_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:23:59
% EndTime: 2019-03-08 18:23:59
% DurationCPUTime: 0.10s
% Computational Cost: add. (90->24), mult. (254->39), div. (0->0), fcn. (191->4), ass. (0->20)
t25 = mrSges(4,1) + mrSges(5,1);
t13 = cos(qJ(3));
t24 = t13 * pkin(2);
t18 = pkin(2) * qJD(3);
t22 = qJD(2) + qJD(3);
t21 = m(5) * pkin(3);
t11 = sin(qJ(3));
t12 = sin(qJ(2));
t14 = cos(qJ(2));
t9 = -t11 * t12 + t13 * t14;
t5 = t22 * t9;
t8 = t11 * t14 + t13 * t12;
t20 = t11 * pkin(2) * t5 + t13 * t8 * t18;
t19 = -mrSges(4,2) - mrSges(5,2);
t17 = qJD(3) * t11 * t9;
t16 = t19 * t13;
t6 = t22 * t8;
t15 = t19 * t5 - t25 * t6;
t10 = pkin(3) + t24;
t1 = [0.4e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * (t8 * t5 - t9 * t6); (-t12 * mrSges(3,1) - t14 * mrSges(3,2)) * qJD(2) + m(4) * ((-t13 * t6 - t17) * pkin(2) + t20) + m(5) * (-pkin(2) * t17 - t10 * t6 + t20) + t15; 0.2e1 * (t16 + ((-t10 + t24) * m(5) - t25) * t11) * t18; -t6 * t21 + t15; (t16 + (-t21 - t25) * t11) * t18; 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
