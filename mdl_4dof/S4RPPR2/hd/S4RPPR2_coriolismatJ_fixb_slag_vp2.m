% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4RPPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:47:26
% EndTime: 2018-11-14 13:47:27
% DurationCPUTime: 0.11s
% Computational Cost: add. (452->24), mult. (625->32), div. (0->0), fcn. (572->4), ass. (0->21)
t20 = sin(pkin(6));
t21 = cos(pkin(6));
t22 = sin(qJ(4));
t23 = cos(qJ(4));
t13 = -t22 * t20 + t23 * t21;
t14 = -t23 * t20 - t22 * t21;
t26 = -t14 * mrSges(5,1) + t13 * mrSges(5,2);
t32 = t26 * qJD(4);
t31 = t26 * qJD(1);
t24 = -pkin(1) - pkin(2);
t25 = -t20 * qJ(2) + t21 * t24;
t15 = -pkin(3) + t25;
t16 = t21 * qJ(2) + t20 * t24;
t10 = t23 * t15 - t22 * t16;
t11 = t22 * t15 + t23 * t16;
t2 = t11 * mrSges(5,1) + t10 * mrSges(5,2);
t30 = t2 * qJD(1);
t29 = t2 * qJD(4);
t5 = m(3) * qJ(2) + t20 * mrSges(4,1) + t21 * mrSges(4,2) + mrSges(3,3) + m(5) * (t10 * t14 + t11 * t13) + m(4) * (t16 * t21 - t25 * t20) + t26;
t28 = t5 * qJD(1);
t1 = [t5 * qJD(2) + t29, t28, 0, -t29 + t30; -t28 + t32, 0, 0, t31 - t32; 0, 0, 0, 0; -qJD(2) * t26 - t30, -t31, 0, 0;];
Cq  = t1;
