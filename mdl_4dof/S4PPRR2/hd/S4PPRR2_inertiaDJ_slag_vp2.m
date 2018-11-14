% Calculate time derivative of joint inertia matrix for
% S4PPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
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
% Datum: 2018-11-14 14:00
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PPRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR2_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR2_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR2_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:59:22
% EndTime: 2018-11-14 13:59:23
% DurationCPUTime: 0.10s
% Computational Cost: add. (87->19), mult. (236->40), div. (0->0), fcn. (232->6), ass. (0->17)
t13 = sin(qJ(4));
t15 = cos(qJ(4));
t11 = sin(pkin(6));
t12 = cos(pkin(6));
t14 = sin(qJ(3));
t16 = cos(qJ(3));
t10 = -t14 * t11 + t16 * t12;
t9 = t16 * t11 + t14 * t12;
t4 = t15 * t10 - t13 * t9;
t6 = t10 * qJD(3);
t7 = t9 * qJD(3);
t2 = qJD(4) * t4 - t13 * t7 + t15 * t6;
t5 = t13 * t10 + t15 * t9;
t3 = -qJD(4) * t5 - t13 * t6 - t15 * t7;
t17 = t3 * mrSges(5,1) - t2 * mrSges(5,2);
t8 = (-mrSges(5,1) * t13 - mrSges(5,2) * t15) * qJD(4) * pkin(3);
t1 = [0.2e1 * m(4) * (-t10 * t7 + t9 * t6) + 0.2e1 * m(5) * (t5 * t2 + t4 * t3); 0; 0; -t7 * mrSges(4,1) - t6 * mrSges(4,2) + m(5) * (t13 * t2 + t15 * t3 + (-t13 * t4 + t15 * t5) * qJD(4)) * pkin(3) + t17; 0; 0.2e1 * t8; t17; 0; t8; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
