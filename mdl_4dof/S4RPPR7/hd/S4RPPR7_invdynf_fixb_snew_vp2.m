% Calculate vector of cutting forces with Newton-Euler
% S4RPPR7
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
% f_new [3x5]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RPPR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR7_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR7_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR7_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:36
% EndTime: 2019-12-31 16:41:37
% DurationCPUTime: 0.30s
% Computational Cost: add. (1582->84), mult. (3356->108), div. (0->0), fcn. (1868->6), ass. (0->49)
t45 = qJD(1) ^ 2;
t42 = sin(qJ(1));
t44 = cos(qJ(1));
t59 = t42 * g(1) - t44 * g(2);
t50 = -t45 * qJ(2) + qJDD(2) - t59;
t66 = -pkin(1) - qJ(3);
t73 = -(2 * qJD(1) * qJD(3)) + qJDD(1) * t66 + t50;
t39 = sin(pkin(6));
t36 = t39 ^ 2;
t40 = cos(pkin(6));
t65 = t40 ^ 2 + t36;
t58 = t65 * mrSges(4,3);
t56 = -t44 * g(1) - t42 * g(2);
t72 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t56;
t71 = -m(2) - m(3);
t70 = pkin(3) * t45;
t69 = mrSges(4,2) * t40;
t68 = mrSges(2,1) - mrSges(3,2);
t67 = -mrSges(2,2) + mrSges(3,3);
t64 = qJDD(1) * t39;
t62 = t39 * g(3) + t73 * t40;
t41 = sin(qJ(4));
t43 = cos(qJ(4));
t52 = -qJDD(1) * mrSges(4,3) - t45 * (mrSges(4,1) * t39 + t69);
t55 = -t39 * t43 - t40 * t41;
t26 = t55 * qJD(1);
t54 = -t39 * t41 + t40 * t43;
t27 = t54 * qJD(1);
t13 = -t26 * mrSges(5,1) + t27 * mrSges(5,2);
t16 = t26 * qJD(4) + qJDD(1) * t54;
t21 = -qJD(4) * mrSges(5,2) + t26 * mrSges(5,3);
t8 = (-pkin(5) * qJDD(1) - t39 * t70) * t40 + t62;
t57 = -t40 * g(3) + t73 * t39;
t9 = -pkin(5) * t64 - t36 * t70 + t57;
t6 = m(5) * (-t41 * t9 + t43 * t8) - t16 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t27 * t13 + qJD(4) * t21;
t15 = -t27 * qJD(4) + qJDD(1) * t55;
t22 = qJD(4) * mrSges(5,1) - t27 * mrSges(5,3);
t7 = m(5) * (t41 * t8 + t43 * t9) + t15 * mrSges(5,3) - qJDD(4) * mrSges(5,2) + t26 * t13 - qJD(4) * t22;
t3 = m(4) * t62 + t40 * t52 + t41 * t7 + t43 * t6;
t4 = m(4) * t57 + t39 * t52 - t41 * t6 + t43 * t7;
t60 = -t39 * t3 + t40 * t4;
t51 = -m(3) * (-qJDD(1) * pkin(1) + t50) - t40 * t3 - t39 * t4;
t48 = qJDD(3) + t72;
t49 = t15 * mrSges(5,1) + t26 * t21 - m(5) * (pkin(3) * t64 + (-pkin(5) * t65 + t66) * t45 + t48) - t27 * t22 - t16 * mrSges(5,2);
t47 = m(4) * (t45 * t66 + t48) + mrSges(4,1) * t64 + qJDD(1) * t69 - t49;
t46 = m(3) * (t45 * pkin(1) - t72) - t47;
t5 = m(2) * t56 + t67 * qJDD(1) + (-t58 - t68) * t45 - t46;
t1 = m(2) * t59 + qJDD(1) * t68 + t45 * t67 + t51;
t2 = [-m(1) * g(1) - t42 * t1 + t44 * t5, t5, -m(3) * g(3) + t60, t4, t7; -m(1) * g(2) + t44 * t1 + t42 * t5, t1, -qJDD(1) * mrSges(3,3) + (-mrSges(3,2) + t58) * t45 + t46, t3, t6; (-m(1) + t71) * g(3) + t60, g(3) * t71 + t60, qJDD(1) * mrSges(3,2) - t45 * mrSges(3,3) - t51, -t45 * t58 + t47, -t49;];
f_new = t2;
