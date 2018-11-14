% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4RPRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:48:28
% EndTime: 2018-11-14 13:48:28
% DurationCPUTime: 0.09s
% Computational Cost: add. (181->26), mult. (374->33), div. (0->0), fcn. (228->4), ass. (0->19)
t14 = cos(qJ(3));
t21 = pkin(1) * sin(pkin(6));
t17 = t14 * t21;
t10 = cos(pkin(6)) * pkin(1) + pkin(2);
t13 = sin(qJ(3));
t20 = t13 * t10;
t9 = t17 + t20;
t7 = qJ(4) + t9;
t24 = m(5) * t7 + mrSges(5,3);
t25 = t24 * qJD(4);
t8 = t14 * t10 - t13 * t21;
t15 = (-mrSges(4,1) - mrSges(5,1)) * t9 + (-mrSges(4,2) + mrSges(5,3)) * t8;
t1 = m(5) * (t7 * t8 + (-pkin(3) - t8) * t9) + t15;
t19 = t1 * qJD(1);
t18 = t24 * qJD(1);
t11 = m(5) * qJ(4) + mrSges(5,3);
t2 = -mrSges(5,3) + 0.2e1 * (t9 / 0.4e1 - t17 / 0.4e1 - t20 / 0.4e1 - qJ(4) / 0.2e1) * m(5);
t16 = t2 * qJD(1) - t11 * qJD(3);
t3 = [t1 * qJD(3) + t25, 0, t19 + (m(5) * (-pkin(3) * t9 + qJ(4) * t8) + t15) * qJD(3) + t25, qJD(3) * t24 + t18; 0, 0, 0, 0; -t2 * qJD(4) - t19, 0, t11 * qJD(4), -t16; t2 * qJD(3) - t18, 0, t16, 0;];
Cq  = t3;
