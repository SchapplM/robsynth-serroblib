% Calculate vector of cutting forces with Newton-Euler
% S4RPRP5
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
% f_new [3x5]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RPRP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP5_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP5_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP5_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:47
% EndTime: 2019-12-31 16:44:48
% DurationCPUTime: 0.38s
% Computational Cost: add. (2296->107), mult. (5449->134), div. (0->0), fcn. (3307->6), ass. (0->56)
t53 = qJD(1) ^ 2;
t48 = cos(pkin(6));
t45 = t48 ^ 2;
t47 = sin(pkin(6));
t73 = t47 ^ 2 + t45;
t83 = t73 * mrSges(3,3);
t49 = sin(qJ(3));
t78 = cos(qJ(3));
t82 = t47 * t49 - t48 * t78;
t50 = sin(qJ(1));
t51 = cos(qJ(1));
t63 = -g(1) * t51 - g(2) * t50;
t37 = -pkin(1) * t53 + qJDD(1) * qJ(2) + t63;
t61 = -mrSges(3,1) * t48 + mrSges(3,2) * t47;
t60 = mrSges(3,3) * qJDD(1) + t53 * t61;
t69 = qJD(1) * qJD(2);
t65 = -g(3) * t48 - 0.2e1 * t47 * t69;
t58 = t78 * t47 + t48 * t49;
t36 = t58 * qJD(1);
t71 = qJD(3) * t36;
t26 = t82 * qJDD(1) + t71;
t30 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t36;
t35 = t82 * qJD(1);
t20 = pkin(3) * t35 - qJ(4) * t36;
t31 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t36;
t52 = qJD(3) ^ 2;
t72 = pkin(5) * qJDD(1);
t79 = pkin(2) * t53;
t15 = (t48 * t79 - t37 - t72) * t47 + t65;
t64 = -g(3) * t47 + (t37 + 0.2e1 * t69) * t48;
t16 = -t45 * t79 + t48 * t72 + t64;
t75 = t49 * t15 + t78 * t16;
t68 = m(5) * (-pkin(3) * t52 + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) - t20 * t35 + t75) + qJD(3) * t31 + qJDD(3) * mrSges(5,3);
t21 = mrSges(5,1) * t35 - mrSges(5,3) * t36;
t74 = -mrSges(4,1) * t35 - mrSges(4,2) * t36 - t21;
t76 = -mrSges(4,3) - mrSges(5,2);
t7 = m(4) * t75 - qJDD(3) * mrSges(4,2) - qJD(3) * t30 + t76 * t26 + t74 * t35 + t68;
t70 = t35 * qJD(3);
t27 = t58 * qJDD(1) - t70;
t29 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t35;
t32 = -mrSges(5,2) * t35 + qJD(3) * mrSges(5,3);
t59 = t78 * t15 - t49 * t16;
t80 = m(5) * (-qJDD(3) * pkin(3) - t52 * qJ(4) + t36 * t20 + qJDD(4) - t59);
t8 = m(4) * t59 - t80 + t74 * t36 + t76 * t27 + (mrSges(4,1) + mrSges(5,1)) * qJDD(3) + (t29 + t32) * qJD(3);
t4 = m(3) * t65 + t49 * t7 + t78 * t8 + (-m(3) * t37 - t60) * t47;
t5 = m(3) * t64 + t60 * t48 - t49 * t8 + t78 * t7;
t81 = t48 * t4 + t47 * t5;
t66 = g(1) * t50 - t51 * g(2);
t62 = qJDD(2) - t66;
t55 = (-pkin(2) * t48 - pkin(1)) * qJDD(1) + (-t73 * pkin(5) - qJ(2)) * t53 + t62;
t57 = t27 * mrSges(5,3) + t36 * t31 - t35 * t32 - t26 * mrSges(5,1) - m(5) * (-0.2e1 * qJD(4) * t36 + (-t27 + t70) * qJ(4) + (t26 + t71) * pkin(3) + t55);
t56 = m(4) * t55 + t26 * mrSges(4,1) + t27 * mrSges(4,2) + t35 * t29 + t36 * t30 - t57;
t54 = m(3) * (-qJDD(1) * pkin(1) - qJ(2) * t53 + t62) + t56;
t6 = (-mrSges(2,2) + t83) * t53 + (mrSges(2,1) - t61) * qJDD(1) + m(2) * t66 - t54;
t1 = m(2) * t63 - t53 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t47 * t4 + t48 * t5;
t2 = [-m(1) * g(1) + t1 * t51 - t50 * t6, t1, t5, t7, -t26 * mrSges(5,2) - t35 * t21 + t68; -m(1) * g(2) + t1 * t50 + t51 * t6, t6, t4, t8, -t57; (-m(1) - m(2)) * g(3) + t81, -m(2) * g(3) + t81, t61 * qJDD(1) - t53 * t83 + t54, t56, -qJDD(3) * mrSges(5,1) + t27 * mrSges(5,2) - qJD(3) * t32 + t36 * t21 + t80;];
f_new = t2;
