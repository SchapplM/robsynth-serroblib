% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
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
% Datum: 2018-11-14 14:09
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4PRPP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP4_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP4_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP4_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP4_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:09:19
% EndTime: 2018-11-14 14:09:19
% DurationCPUTime: 0.07s
% Computational Cost: add. (62->17), mult. (158->26), div. (0->0), fcn. (145->4), ass. (0->11)
t11 = sin(pkin(5));
t12 = cos(pkin(5));
t15 = sin(qJ(2));
t16 = cos(qJ(2));
t4 = -t11 * t16 - t12 * t15;
t17 = m(5) * t4;
t6 = t11 * pkin(2) + qJ(4);
t5 = m(5) * t6 + mrSges(5,3);
t13 = t5 * qJD(2);
t3 = t11 * t15 - t12 * t16;
t1 = [0 (-t16 * mrSges(3,2) - t15 * mrSges(3,1) + t3 * mrSges(4,2) + t4 * mrSges(4,1) + m(4) * (-t11 * t3 + t12 * t4) * pkin(2) + t4 * mrSges(5,1) - t3 * mrSges(5,3) + m(5) * (-t6 * t3 - (-t12 * pkin(2) - pkin(3)) * t4)) * qJD(2) - qJD(4) * t17, 0, -qJD(2) * t17; 0, t5 * qJD(4), 0, t13; 0, 0, 0, 0; 0, -t13, 0, 0;];
Cq  = t1;
