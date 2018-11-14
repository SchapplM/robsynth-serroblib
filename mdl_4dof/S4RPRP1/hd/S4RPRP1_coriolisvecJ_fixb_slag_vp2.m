% Calculate vector of centrifugal and coriolis load on the joints for
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
% tauc [4x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4RPRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP1_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:48:27
% EndTime: 2018-11-14 13:48:28
% DurationCPUTime: 0.13s
% Computational Cost: add. (159->38), mult. (409->52), div. (0->0), fcn. (180->4), ass. (0->26)
t32 = mrSges(4,1) + mrSges(5,1);
t18 = cos(pkin(6)) * pkin(1) + pkin(2);
t22 = sin(qJ(3));
t23 = cos(qJ(3));
t33 = pkin(1) * sin(pkin(6));
t31 = t22 * t18 + t23 * t33;
t15 = t18 * qJD(1);
t27 = qJD(1) * t33;
t26 = t22 * t27;
t8 = t23 * t15 - t26;
t30 = qJD(4) - t8;
t29 = qJD(3) * t23;
t28 = t22 * t33;
t11 = -qJD(3) * t28 + t18 * t29;
t25 = t23 * t18 - t28;
t19 = qJD(1) + qJD(3);
t6 = -qJD(3) * t26 + t15 * t29;
t2 = t19 * qJD(4) + t6;
t9 = t22 * t15 + t23 * t27;
t7 = t9 * qJD(3);
t24 = -t6 * mrSges(4,2) + t2 * mrSges(5,3) - t32 * t7;
t12 = t31 * qJD(3);
t10 = qJD(4) + t11;
t4 = t19 * qJ(4) + t9;
t3 = -t19 * pkin(3) + t30;
t1 = [m(4) * (t9 * t11 - t8 * t12 - t7 * t25 + t6 * t31) + m(5) * (t2 * (qJ(4) + t31) + t4 * t10 + t7 * (-pkin(3) - t25) + t3 * t12) + (-t11 * mrSges(4,2) + t10 * mrSges(5,3) - t32 * t12) * t19 + t24; 0; (t8 * mrSges(4,2) + t30 * mrSges(5,3) + t32 * t9) * t19 + t24 + (-t7 * pkin(3) + t2 * qJ(4) - t3 * t9 + t30 * t4) * m(5); -t19 ^ 2 * mrSges(5,3) + (-t4 * t19 + t7) * m(5);];
tauc  = t1(:);
