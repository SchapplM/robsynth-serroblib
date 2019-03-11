% Calculate time derivative of joint inertia matrix for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [3x3]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S3RRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_inertiaDJ_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_inertiaDJ_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_inertiaDJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_inertiaDJ_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRR1_inertiaDJ_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RRR1_inertiaDJ_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:07:57
% EndTime: 2019-03-08 18:07:57
% DurationCPUTime: 0.09s
% Computational Cost: add. (65->21), mult. (200->43), div. (0->0), fcn. (110->4), ass. (0->17)
t11 = cos(qJ(2));
t10 = cos(qJ(3));
t14 = qJD(3) * t10;
t8 = sin(qJ(3));
t15 = qJD(3) * t8;
t9 = sin(qJ(2));
t18 = t8 * t9;
t7 = t11 * pkin(1) + pkin(2);
t2 = t7 * t14 + (-t9 * t15 + (t10 * t11 - t18) * qJD(2)) * pkin(1);
t17 = t10 * t9;
t3 = -t7 * t15 + (-t9 * t14 + (-t11 * t8 - t17) * qJD(2)) * pkin(1);
t13 = t3 * mrSges(4,1) - t2 * mrSges(4,2);
t19 = (-mrSges(3,1) * t9 - mrSges(3,2) * t11) * qJD(2) * pkin(1) + t13;
t6 = (-t8 * mrSges(4,1) - t10 * mrSges(4,2)) * qJD(3) * pkin(2);
t5 = pkin(1) * t17 + t8 * t7;
t4 = -pkin(1) * t18 + t10 * t7;
t1 = [0.2e1 * m(4) * (t5 * t2 + t4 * t3) + 0.2e1 * t19; (m(4) * (t10 * t3 + t5 * t14 - t4 * t15 + t2 * t8) - mrSges(4,2) * t14 - mrSges(4,1) * t15) * pkin(2) + t19; 0.2e1 * t6; t13; t6; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t1(1) t1(2) t1(4); t1(2) t1(3) t1(5); t1(4) t1(5) t1(6);];
Mq  = res;
