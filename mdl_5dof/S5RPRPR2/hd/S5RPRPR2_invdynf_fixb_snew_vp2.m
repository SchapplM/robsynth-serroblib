% Calculate vector of cutting forces with Newton-Euler
% S5RPRPR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:25
% EndTime: 2019-12-05 17:49:26
% DurationCPUTime: 0.70s
% Computational Cost: add. (9552->100), mult. (13789->133), div. (0->0), fcn. (7969->10), ass. (0->63)
t53 = qJD(1) + qJD(3);
t49 = t53 ^ 2;
t57 = cos(pkin(9));
t52 = t57 ^ 2;
t55 = sin(pkin(9));
t80 = t55 ^ 2 + t52;
t87 = t80 * mrSges(5,3);
t86 = pkin(4) * t57;
t50 = qJDD(1) + qJDD(3);
t85 = pkin(7) * t50;
t61 = sin(qJ(1));
t64 = cos(qJ(1));
t81 = t64 * g(2) + t61 * g(3);
t38 = qJDD(1) * pkin(1) + t81;
t65 = qJD(1) ^ 2;
t76 = t61 * g(2) - t64 * g(3);
t39 = -t65 * pkin(1) + t76;
t56 = sin(pkin(8));
t58 = cos(pkin(8));
t73 = t58 * t38 - t56 * t39;
t28 = qJDD(1) * pkin(2) + t73;
t83 = t56 * t38 + t58 * t39;
t29 = -t65 * pkin(2) + t83;
t60 = sin(qJ(3));
t63 = cos(qJ(3));
t84 = t60 * t28 + t63 * t29;
t54 = -g(1) + qJDD(2);
t79 = qJD(4) * t53;
t82 = t57 * t54 - 0.2e1 * t55 * t79;
t19 = -t49 * pkin(3) + t50 * qJ(4) + t84;
t13 = (t49 * t86 - t19 - t85) * t55 + t82;
t77 = t55 * t54 + (t19 + 0.2e1 * t79) * t57;
t14 = -t52 * t49 * pkin(4) + t57 * t85 + t77;
t59 = sin(qJ(5));
t62 = cos(qJ(5));
t68 = -t55 * t59 + t57 * t62;
t32 = t68 * t53;
t69 = t55 * t62 + t57 * t59;
t33 = t69 * t53;
t22 = -t32 * mrSges(6,1) + t33 * mrSges(6,2);
t24 = t32 * qJD(5) + t69 * t50;
t30 = -qJD(5) * mrSges(6,2) + t32 * mrSges(6,3);
t11 = m(6) * (t62 * t13 - t59 * t14) - t24 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t33 * t22 + qJD(5) * t30;
t23 = -t33 * qJD(5) + t68 * t50;
t31 = qJD(5) * mrSges(6,1) - t33 * mrSges(6,3);
t12 = m(6) * (t59 * t13 + t62 * t14) + t23 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t32 * t22 - qJD(5) * t31;
t71 = -t57 * mrSges(5,1) + t55 * mrSges(5,2);
t70 = t50 * mrSges(5,3) + t49 * t71;
t8 = m(5) * t82 + t59 * t12 + t62 * t11 + (-m(5) * t19 - t70) * t55;
t9 = m(5) * t77 - t59 * t11 + t62 * t12 + t70 * t57;
t78 = m(4) * t54 + t55 * t9 + t57 * t8;
t75 = m(3) * t54 + t78;
t74 = t63 * t28 - t60 * t29;
t72 = qJDD(4) - t74;
t67 = t23 * mrSges(6,1) + t32 * t30 - m(6) * ((-pkin(3) - t86) * t50 + (-t80 * pkin(7) - qJ(4)) * t49 + t72) - t33 * t31 - t24 * mrSges(6,2);
t66 = m(5) * (-t50 * pkin(3) - t49 * qJ(4) + t72) - t67;
t10 = m(4) * t74 + (mrSges(4,1) - t71) * t50 + (-mrSges(4,2) + t87) * t49 - t66;
t5 = m(4) * t84 - t49 * mrSges(4,1) - t50 * mrSges(4,2) - t55 * t8 + t57 * t9;
t4 = m(3) * t83 - t65 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t60 * t10 + t63 * t5;
t3 = m(3) * t73 + qJDD(1) * mrSges(3,1) - t65 * mrSges(3,2) + t63 * t10 + t60 * t5;
t2 = m(2) * t81 + qJDD(1) * mrSges(2,1) - t65 * mrSges(2,2) + t58 * t3 + t56 * t4;
t1 = m(2) * t76 - t65 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t56 * t3 + t58 * t4;
t6 = [(-m(1) - m(2)) * g(1) + t75, t1, t4, t5, t9, t12; -m(1) * g(2) - t61 * t1 - t64 * t2, t2, t3, t10, t8, t11; -m(1) * g(3) + t64 * t1 - t61 * t2, -m(2) * g(1) + t75, t75, t78, -t49 * t87 + t71 * t50 + t66, -t67;];
f_new = t6;
