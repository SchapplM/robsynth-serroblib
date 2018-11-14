% Calculate time derivative of joint inertia matrix for
% S4PRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:02
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PRPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR2_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR2_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR2_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:02:11
% EndTime: 2018-11-14 14:02:11
% DurationCPUTime: 0.08s
% Computational Cost: add. (119->28), mult. (312->55), div. (0->0), fcn. (284->6), ass. (0->23)
t15 = sin(pkin(6));
t23 = pkin(2) * t15;
t16 = cos(pkin(6));
t18 = sin(qJ(2));
t20 = cos(qJ(2));
t12 = t15 * t20 + t16 * t18;
t10 = t12 * qJD(2);
t13 = -t15 * t18 + t16 * t20;
t11 = t13 * qJD(2);
t17 = sin(qJ(4));
t19 = cos(qJ(4));
t4 = -t17 * t12 + t19 * t13;
t2 = t4 * qJD(4) - t17 * t10 + t19 * t11;
t5 = t19 * t12 + t17 * t13;
t3 = -t5 * qJD(4) - t19 * t10 - t17 * t11;
t22 = t3 * mrSges(5,1) - t2 * mrSges(5,2);
t14 = t16 * pkin(2) + pkin(3);
t8 = t19 * t14 - t17 * t23;
t6 = t8 * qJD(4);
t9 = t17 * t14 + t19 * t23;
t7 = t9 * qJD(4);
t21 = -t7 * mrSges(5,1) - t6 * mrSges(5,2);
t1 = [0.2e1 * m(4) * (-t13 * t10 + t12 * t11) + 0.2e1 * m(5) * (t5 * t2 + t4 * t3); -t10 * mrSges(4,1) - t11 * mrSges(4,2) + (-t18 * mrSges(3,1) - t20 * mrSges(3,2)) * qJD(2) + m(5) * (t9 * t2 + t8 * t3 - t7 * t4 + t6 * t5) + m(4) * (-t10 * t16 + t11 * t15) * pkin(2) + t22; 0.2e1 * m(5) * (t9 * t6 - t8 * t7) + 0.2e1 * t21; 0; 0; 0; t22; t21; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
