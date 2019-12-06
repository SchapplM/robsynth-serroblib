% Calculate vector of cutting forces with Newton-Euler
% S5RRRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_invdynf_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:52:58
% EndTime: 2019-12-05 18:53:01
% DurationCPUTime: 0.87s
% Computational Cost: add. (12258->121), mult. (16543->163), div. (0->0), fcn. (11435->10), ass. (0->70)
t54 = qJD(1) + qJD(2);
t56 = sin(qJ(4));
t57 = sin(qJ(3));
t61 = cos(qJ(4));
t62 = cos(qJ(3));
t40 = (t56 * t57 - t61 * t62) * t54;
t52 = qJDD(1) + qJDD(2);
t76 = qJD(3) * t54;
t43 = t57 * t52 + t62 * t76;
t75 = t57 * t76;
t44 = t62 * t52 - t75;
t78 = t54 * t57;
t47 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t78;
t77 = t54 * t62;
t48 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t77;
t50 = t54 ^ 2;
t59 = sin(qJ(1));
t64 = cos(qJ(1));
t74 = t59 * g(1) - t64 * g(2);
t45 = qJDD(1) * pkin(1) + t74;
t65 = qJD(1) ^ 2;
t72 = -t64 * g(1) - t59 * g(2);
t46 = -t65 * pkin(1) + t72;
t58 = sin(qJ(2));
t63 = cos(qJ(2));
t34 = t58 * t45 + t63 * t46;
t71 = -t62 * g(3) - t57 * t34;
t28 = (t50 * t57 * t62 + qJDD(3)) * pkin(2) + t71;
t73 = -t57 * g(3) + t62 * t34;
t29 = (-t50 * t62 ^ 2 - qJD(3) ^ 2) * pkin(2) + t73;
t16 = t56 * t28 + t61 * t29;
t23 = -t40 * qJD(4) + t61 * t43 + t56 * t44;
t41 = (t56 * t62 + t57 * t61) * t54;
t53 = qJD(3) + qJD(4);
t55 = sin(qJ(5));
t60 = cos(qJ(5));
t35 = -t55 * t41 + t60 * t53;
t51 = qJDD(3) + qJDD(4);
t18 = t35 * qJD(5) + t60 * t23 + t55 * t51;
t36 = t60 * t41 + t55 * t53;
t19 = -t35 * mrSges(6,1) + t36 * mrSges(6,2);
t22 = -t41 * qJD(4) - t56 * t43 + t61 * t44;
t21 = qJDD(5) - t22;
t33 = -t63 * t45 + t58 * t46;
t25 = (-t44 + t75) * pkin(2) + t33;
t39 = qJD(5) + t40;
t26 = -t39 * mrSges(6,2) + t35 * mrSges(6,3);
t13 = m(6) * (-t55 * t16 + t60 * t25) - t18 * mrSges(6,3) + t21 * mrSges(6,1) - t36 * t19 + t39 * t26;
t17 = -t36 * qJD(5) - t55 * t23 + t60 * t51;
t27 = t39 * mrSges(6,1) - t36 * mrSges(6,3);
t14 = m(6) * (t60 * t16 + t55 * t25) + t17 * mrSges(6,3) - t21 * mrSges(6,2) + t35 * t19 - t39 * t27;
t37 = -t53 * mrSges(5,2) - t40 * mrSges(5,3);
t38 = t53 * mrSges(5,1) - t41 * mrSges(5,3);
t68 = -m(5) * t25 + t22 * mrSges(5,1) - t23 * mrSges(5,2) - t60 * t13 - t55 * t14 - t40 * t37 - t41 * t38;
t81 = t44 * mrSges(4,1) - t43 * mrSges(4,2) - (t47 * t57 - t48 * t62) * t54 + t68;
t80 = -m(2) - m(3);
t15 = -t61 * t28 + t56 * t29;
t31 = t40 * mrSges(5,1) + t41 * mrSges(5,2);
t67 = t17 * mrSges(6,1) - t18 * mrSges(6,2) + t35 * t26 - t36 * t27;
t10 = t51 * mrSges(5,1) - t23 * mrSges(5,3) - t41 * t31 + t53 * t37 + (-m(5) - m(6)) * t15 + t67;
t42 = (-mrSges(4,1) * t62 + mrSges(4,2) * t57) * t54;
t9 = m(5) * t16 - t51 * mrSges(5,2) + t22 * mrSges(5,3) - t55 * t13 + t60 * t14 - t40 * t31 - t53 * t38;
t6 = m(4) * t71 + qJDD(3) * mrSges(4,1) - t43 * mrSges(4,3) + qJD(3) * t48 + t61 * t10 - t42 * t78 + t56 * t9;
t7 = m(4) * t73 - qJDD(3) * mrSges(4,2) + t44 * mrSges(4,3) - qJD(3) * t47 - t56 * t10 + t42 * t77 + t61 * t9;
t79 = t57 * t7 + t62 * t6;
t8 = t52 * mrSges(3,1) - t50 * mrSges(3,2) + (-m(3) - m(4)) * t33 + t81;
t3 = m(3) * t34 - t50 * mrSges(3,1) - t52 * mrSges(3,2) - t57 * t6 + t62 * t7;
t2 = m(2) * t72 - t65 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t63 * t3 - t58 * t8;
t1 = m(2) * t74 + qJDD(1) * mrSges(2,1) - t65 * mrSges(2,2) + t58 * t3 + t63 * t8;
t4 = [-m(1) * g(1) - t59 * t1 + t64 * t2, t2, t3, t7, t9, t14; -m(1) * g(2) + t64 * t1 + t59 * t2, t1, t8, t6, t10, t13; (-m(1) + t80) * g(3) + t79, t80 * g(3) + t79, -m(3) * g(3) + t79, m(4) * t33 - t81, -t68, m(6) * t15 - t67;];
f_new = t4;
