% Calculate time derivative of joint inertia matrix for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-03-08 18:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR2_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR2_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:28:26
% EndTime: 2019-03-08 18:28:27
% DurationCPUTime: 0.08s
% Computational Cost: add. (121->25), mult. (215->44), div. (0->0), fcn. (160->4), ass. (0->19)
t12 = sin(pkin(6));
t13 = cos(pkin(6));
t16 = -pkin(1) - pkin(2);
t19 = -t12 * qJ(2) + t13 * t16;
t10 = t13 * qJ(2) + t12 * t16;
t14 = sin(qJ(4));
t15 = cos(qJ(4));
t9 = -pkin(3) + t19;
t3 = -t14 * t10 + t15 * t9;
t7 = -t14 * t12 + t15 * t13;
t1 = t7 * qJD(2) + t3 * qJD(4);
t4 = t15 * t10 + t14 * t9;
t8 = t15 * t12 + t14 * t13;
t2 = -t8 * qJD(2) - t4 * qJD(4);
t18 = t2 * mrSges(5,1) - t1 * mrSges(5,2);
t5 = t7 * qJD(4);
t6 = t8 * qJD(4);
t17 = -t6 * mrSges(5,1) - t5 * mrSges(5,2);
t11 = [0.2e1 * m(5) * (t4 * t1 + t3 * t2) - 0.2e1 * t18 + 0.2e1 * (t13 * mrSges(4,2) + t12 * mrSges(4,1) + m(3) * qJ(2) + mrSges(3,3) + m(4) * (t10 * t13 - t19 * t12)) * qJD(2); m(5) * (t8 * t1 + t7 * t2 - t6 * t3 + t5 * t4) - t17; 0.2e1 * m(5) * (t8 * t5 - t7 * t6); 0; 0; 0; t18; t17; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t11(1) t11(2) t11(4) t11(7); t11(2) t11(3) t11(5) t11(8); t11(4) t11(5) t11(6) t11(9); t11(7) t11(8) t11(9) t11(10);];
Mq  = res;
