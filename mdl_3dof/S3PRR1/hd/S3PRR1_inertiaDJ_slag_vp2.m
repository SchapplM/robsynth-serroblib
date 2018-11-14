% Calculate time derivative of joint inertia matrix for
% S3PRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S3PRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_inertiaDJ_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_inertiaDJ_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_inertiaDJ_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR1_inertiaDJ_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PRR1_inertiaDJ_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PRR1_inertiaDJ_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:11:57
% EndTime: 2018-11-14 10:11:57
% DurationCPUTime: 0.06s
% Computational Cost: add. (39->12), mult. (116->27), div. (0->0), fcn. (90->4), ass. (0->12)
t12 = qJD(2) + qJD(3);
t10 = cos(qJ(2));
t7 = sin(qJ(3));
t8 = sin(qJ(2));
t9 = cos(qJ(3));
t6 = t9 * t10 - t7 * t8;
t2 = t12 * t6;
t5 = t7 * t10 + t9 * t8;
t3 = t12 * t5;
t11 = -t3 * mrSges(4,1) - t2 * mrSges(4,2);
t4 = (-mrSges(4,1) * t7 - mrSges(4,2) * t9) * qJD(3) * pkin(2);
t1 = [0.2e1 * m(4) * (t5 * t2 - t6 * t3); (-t8 * mrSges(3,1) - t10 * mrSges(3,2)) * qJD(2) + m(4) * (t2 * t7 - t3 * t9 + (t5 * t9 - t6 * t7) * qJD(3)) * pkin(2) + t11; 0.2e1 * t4; t11; t4; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t1(1) t1(2) t1(4); t1(2) t1(3) t1(5); t1(4) t1(5) t1(6);];
Mq  = res;
