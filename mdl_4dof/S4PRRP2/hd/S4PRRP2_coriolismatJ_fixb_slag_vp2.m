% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:03
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4PRRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP2_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:03:19
% EndTime: 2018-11-14 14:03:19
% DurationCPUTime: 0.11s
% Computational Cost: add. (201->20), mult. (503->28), div. (0->0), fcn. (455->4), ass. (0->19)
t17 = sin(qJ(3));
t18 = sin(qJ(2));
t19 = cos(qJ(3));
t20 = cos(qJ(2));
t15 = -t17 * t20 - t19 * t18;
t16 = -t17 * t18 + t19 * t20;
t29 = mrSges(4,2) + mrSges(5,2);
t30 = mrSges(4,1) + mrSges(5,1);
t31 = t30 * t15 - t29 * t16;
t27 = m(5) * pkin(3);
t28 = t27 + t30;
t26 = pkin(2) * t17;
t25 = t19 * pkin(2);
t1 = t15 * t27 + t31;
t22 = t1 * qJD(3);
t7 = t29 * t25 + t28 * t26;
t21 = t7 * qJD(2);
t11 = t16 * t26;
t2 = [0, t22 + (-t18 * mrSges(3,1) - t20 * mrSges(3,2) + m(4) * (t15 * t25 + t11) + m(5) * (t11 + (pkin(3) + t25) * t15) + t31) * qJD(2), t1 * qJD(2) + t22, 0; 0, -t7 * qJD(3) (-t28 * t17 - t19 * t29) * qJD(3) * pkin(2) - t21, 0; 0, t21, 0, 0; 0, 0, 0, 0;];
Cq  = t2;
