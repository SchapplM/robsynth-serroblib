% Calculate time derivative of joint inertia matrix for
% S4RPRP2
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:50
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4RPRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_inertiaDJ_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP2_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP2_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:49:23
% EndTime: 2018-11-14 13:49:23
% DurationCPUTime: 0.09s
% Computational Cost: add. (97->25), mult. (175->36), div. (0->0), fcn. (84->2), ass. (0->18)
t18 = mrSges(4,2) + mrSges(5,2);
t12 = cos(qJ(3));
t11 = sin(qJ(3));
t13 = -pkin(1) - pkin(2);
t8 = -t11 * qJ(2) + t12 * t13;
t4 = t12 * qJD(2) + qJD(3) * t8;
t22 = t18 * t4;
t21 = t18 * t12;
t9 = t12 * qJ(2) + t11 * t13;
t20 = t9 * qJD(3);
t19 = mrSges(4,1) + mrSges(5,1);
t16 = qJD(3) * t11;
t5 = -t11 * qJD(2) - t20;
t15 = t11 * t4 + (t20 + t5) * t12;
t14 = m(5) * pkin(3) + t19;
t7 = -pkin(3) + t8;
t1 = t9 * t4;
t2 = [0.2e1 * m(4) * (t8 * t5 + t1) + 0.2e1 * m(5) * (t7 * t5 + t1) - 0.2e1 * t19 * t5 + 0.2e1 * t22 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * qJD(2); (t19 * t11 + t21) * qJD(3) + m(4) * (-t8 * t16 + t15) + m(5) * (-t7 * t16 + t15); 0; t14 * t5 - t22; (-t14 * t11 - t21) * qJD(3); 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1) t2(2) t2(4) t2(7); t2(2) t2(3) t2(5) t2(8); t2(4) t2(5) t2(6) t2(9); t2(7) t2(8) t2(9) t2(10);];
Mq  = res;
