% Calculate vector of cutting forces with Newton-Euler
% S4RPRR9
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% f_new [3x5]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RPRR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR9_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR9_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR9_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:09
% EndTime: 2019-12-31 16:56:10
% DurationCPUTime: 0.30s
% Computational Cost: add. (1794->97), mult. (3370->123), div. (0->0), fcn. (1717->6), ass. (0->54)
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t55 = -t45 * g(1) - t42 * g(2);
t70 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t55;
t69 = -m(2) - m(3);
t68 = -pkin(1) - pkin(5);
t41 = sin(qJ(3));
t67 = t41 * g(3);
t66 = (mrSges(2,1) - mrSges(3,2));
t65 = -mrSges(2,2) + mrSges(3,3);
t44 = cos(qJ(3));
t64 = qJD(1) * t44;
t63 = t41 * qJD(1);
t62 = qJD(1) * qJD(3);
t30 = (mrSges(4,1) * t41 + mrSges(4,2) * t44) * qJD(1);
t56 = t44 * t62;
t32 = -t41 * qJDD(1) - t56;
t35 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t64;
t40 = sin(qJ(4));
t43 = cos(qJ(4));
t47 = qJD(1) ^ 2;
t59 = t42 * g(1) - t45 * g(2);
t52 = -t47 * qJ(2) + qJDD(2) - t59;
t21 = t68 * qJDD(1) + t52;
t58 = -t44 * g(3) + t41 * t21;
t57 = t41 * t62;
t33 = t44 * qJDD(1) - t57;
t49 = t68 * t47 - t70;
t10 = (-t33 + t57) * pkin(6) + (-t32 + t56) * pkin(3) + t49;
t31 = (pkin(3) * t41 - pkin(6) * t44) * qJD(1);
t46 = qJD(3) ^ 2;
t12 = -t46 * pkin(3) + qJDD(3) * pkin(6) - t31 * t63 + t58;
t28 = t43 * qJD(3) - t40 * t64;
t14 = t28 * qJD(4) + t40 * qJDD(3) + t43 * t33;
t29 = t40 * qJD(3) + t43 * t64;
t15 = -t28 * mrSges(5,1) + t29 * mrSges(5,2);
t36 = qJD(4) + t63;
t16 = -t36 * mrSges(5,2) + t28 * mrSges(5,3);
t27 = qJDD(4) - t32;
t8 = m(5) * (t43 * t10 - t40 * t12) - t14 * mrSges(5,3) + t27 * mrSges(5,1) - t29 * t15 + t36 * t16;
t13 = -t29 * qJD(4) + t43 * qJDD(3) - t40 * t33;
t17 = t36 * mrSges(5,1) - t29 * mrSges(5,3);
t9 = m(5) * (t40 * t10 + t43 * t12) + t13 * mrSges(5,3) - t27 * mrSges(5,2) + t28 * t15 - t36 * t17;
t4 = m(4) * t58 - qJDD(3) * mrSges(4,2) + t32 * mrSges(4,3) - qJD(3) * t35 - t30 * t63 - t40 * t8 + t43 * t9;
t34 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t63;
t48 = m(5) * (-qJDD(3) * pkin(3) - t46 * pkin(6) - t67 + (qJD(1) * t31 - t21) * t44) - t13 * mrSges(5,1) + t14 * mrSges(5,2) - t28 * t16 + t29 * t17;
t5 = m(4) * (t44 * t21 + t67) - t33 * mrSges(4,3) + qJDD(3) * mrSges(4,1) - t30 * t64 + qJD(3) * t34 - t48;
t60 = t44 * t4 - t41 * t5;
t53 = -m(3) * (-qJDD(1) * pkin(1) + t52) - t41 * t4 - t44 * t5;
t51 = m(4) * t49 - t32 * mrSges(4,1) + t33 * mrSges(4,2) + t34 * t63 + t35 * t64 + t40 * t9 + t43 * t8;
t50 = -m(3) * (t47 * pkin(1) + t70) + t51;
t2 = m(2) * t55 + t65 * qJDD(1) - (t66 * t47) + t50;
t1 = m(2) * t59 + t66 * qJDD(1) + t65 * t47 + t53;
t3 = [-m(1) * g(1) - t42 * t1 + t45 * t2, t2, -m(3) * g(3) + t60, t4, t9; -m(1) * g(2) + t45 * t1 + t42 * t2, t1, -(t47 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t50, t5, t8; (-m(1) + t69) * g(3) + t60, t69 * g(3) + t60, qJDD(1) * mrSges(3,2) - t47 * mrSges(3,3) - t53, t51, t48;];
f_new = t3;
