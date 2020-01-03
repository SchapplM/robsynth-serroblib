% Calculate vector of cutting forces with Newton-Euler
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPRR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR9_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR9_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR9_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:21
% EndTime: 2019-12-31 18:02:22
% DurationCPUTime: 0.49s
% Computational Cost: add. (4694->113), mult. (8234->143), div. (0->0), fcn. (3747->8), ass. (0->62)
t54 = qJD(1) ^ 2;
t49 = sin(qJ(1));
t52 = cos(qJ(1));
t63 = -t52 * g(1) - t49 * g(2);
t60 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t63;
t75 = -pkin(1) - pkin(2);
t24 = t75 * t54 + t60;
t67 = t49 * g(1) - t52 * g(2);
t57 = -t54 * qJ(2) + qJDD(2) - t67;
t26 = t75 * qJDD(1) + t57;
t45 = sin(pkin(8));
t46 = cos(pkin(8));
t64 = -t45 * t24 + t46 * t26;
t14 = qJDD(1) * pkin(3) - t54 * pkin(6) - t64;
t48 = sin(qJ(4));
t51 = cos(qJ(4));
t68 = qJD(1) * qJD(4);
t65 = t51 * t68;
t33 = -t48 * qJDD(1) - t65;
t66 = t48 * t68;
t34 = -t51 * qJDD(1) + t66;
t70 = qJD(1) * t48;
t35 = (qJD(4) * mrSges(5,1)) + mrSges(5,3) * t70;
t69 = t51 * qJD(1);
t36 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t69;
t47 = sin(qJ(5));
t50 = cos(qJ(5));
t10 = (-t33 + t65) * pkin(7) + (-t34 - t66) * pkin(4) + t14;
t32 = (pkin(4) * t51 + pkin(7) * t48) * qJD(1);
t53 = qJD(4) ^ 2;
t71 = t46 * t24 + t45 * t26;
t15 = -t54 * pkin(3) - qJDD(1) * pkin(6) + t71;
t43 = g(3) + qJDD(3);
t72 = t51 * t15 + t48 * t43;
t12 = -t53 * pkin(4) + qJDD(4) * pkin(7) - t32 * t69 + t72;
t29 = t50 * qJD(4) + t47 * t70;
t17 = t29 * qJD(5) + t47 * qJDD(4) + t50 * t33;
t30 = t47 * qJD(4) - t50 * t70;
t18 = -t29 * mrSges(6,1) + t30 * mrSges(6,2);
t37 = qJD(5) + t69;
t22 = -t37 * mrSges(6,2) + t29 * mrSges(6,3);
t28 = qJDD(5) - t34;
t8 = m(6) * (t50 * t10 - t47 * t12) - t17 * mrSges(6,3) + t28 * mrSges(6,1) - t30 * t18 + t37 * t22;
t16 = -t30 * qJD(5) + t50 * qJDD(4) - t47 * t33;
t23 = t37 * mrSges(6,1) - t30 * mrSges(6,3);
t9 = m(6) * (t47 * t10 + t50 * t12) + t16 * mrSges(6,3) - t28 * mrSges(6,2) + t29 * t18 - t37 * t23;
t77 = m(5) * t14 - t34 * mrSges(5,1) + t33 * mrSges(5,2) - (t48 * t35 - t51 * t36) * qJD(1) + t47 * t9 + t50 * t8;
t76 = -m(2) - m(3);
t74 = t51 * t43;
t73 = mrSges(2,1) + mrSges(3,1);
t31 = (mrSges(5,1) * t51 - mrSges(5,2) * t48) * qJD(1);
t6 = m(5) * t72 - qJDD(4) * mrSges(5,2) + t34 * mrSges(5,3) - qJD(4) * t35 - t31 * t69 - t47 * t8 + t50 * t9;
t55 = m(6) * (-qJDD(4) * pkin(4) - t53 * pkin(7) - t74 + (-qJD(1) * t32 + t15) * t48) - t16 * mrSges(6,1) + t17 * mrSges(6,2) - t29 * t22 + t30 * t23;
t7 = m(5) * (-t48 * t15 + t74) - t33 * mrSges(5,3) + qJDD(4) * mrSges(5,1) + t31 * t70 + qJD(4) * t36 - t55;
t4 = m(4) * t71 - t54 * mrSges(4,1) + qJDD(1) * mrSges(4,2) - t48 * t7 + t51 * t6;
t5 = m(4) * t64 - qJDD(1) * mrSges(4,1) - t54 * mrSges(4,2) - t77;
t61 = t46 * t4 - t45 * t5 + m(3) * (-t54 * pkin(1) + t60) + qJDD(1) * mrSges(3,3);
t59 = -m(3) * (-qJDD(1) * pkin(1) + t57) - t45 * t4 - t46 * t5;
t58 = m(4) * t43 + t48 * t6 + t51 * t7;
t2 = m(2) * t67 + (-mrSges(2,2) + mrSges(3,3)) * t54 + t73 * qJDD(1) + t59;
t1 = m(2) * t63 - qJDD(1) * mrSges(2,2) - t73 * t54 + t61;
t3 = [-m(1) * g(1) + t52 * t1 - t49 * t2, t1, -t54 * mrSges(3,1) + t61, t4, t6, t9; -m(1) * g(2) + t49 * t1 + t52 * t2, t2, -m(3) * g(3) - t58, t5, t7, t8; (-m(1) + t76) * g(3) - t58, t76 * g(3) - t58, -qJDD(1) * mrSges(3,1) - t54 * mrSges(3,3) - t59, t58, t77, t55;];
f_new = t3;
