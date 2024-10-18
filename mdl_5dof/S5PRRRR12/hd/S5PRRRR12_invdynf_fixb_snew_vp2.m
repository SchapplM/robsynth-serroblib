% Calculate vector of cutting forces with Newton-Euler
% S5PRRRR12
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRR12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR12_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR12_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR12_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:37
% EndTime: 2024-09-28 18:07:40
% DurationCPUTime: 1.36s
% Computational Cost: add. (26125->89), mult. (35514->134), div. (0->0), fcn. (27110->14), ass. (0->69)
t54 = qJD(2) + qJD(3);
t50 = qJD(4) + t54;
t48 = t50 ^ 2;
t53 = qJDD(2) + qJDD(3);
t49 = qJDD(4) + t53;
t56 = sin(pkin(11));
t59 = cos(pkin(11));
t45 = -g(1) * t59 - g(2) * t56;
t55 = -g(3) + qJDD(1);
t65 = sin(qJ(2));
t69 = cos(qJ(2));
t58 = sin(pkin(5));
t86 = t58 * t69;
t44 = g(1) * t56 - g(2) * t59;
t61 = cos(pkin(5));
t89 = t44 * t61;
t73 = -t45 * t65 + t55 * t86 + t69 * t89;
t29 = qJDD(2) * pkin(2) + t73;
t70 = qJD(2) ^ 2;
t87 = t58 * t65;
t78 = t69 * t45 + t55 * t87 + t65 * t89;
t30 = -pkin(2) * t70 + t78;
t64 = sin(qJ(3));
t68 = cos(qJ(3));
t75 = t68 * t29 - t30 * t64;
t24 = pkin(3) * t53 + t75;
t52 = t54 ^ 2;
t82 = t64 * t29 + t68 * t30;
t25 = -pkin(3) * t52 + t82;
t63 = sin(qJ(4));
t67 = cos(qJ(4));
t83 = t63 * t24 + t67 * t25;
t57 = sin(pkin(6));
t90 = pkin(10) * t57;
t20 = -pkin(4) * t48 + t49 * t90 + t83;
t60 = cos(pkin(6));
t43 = t50 * t60 + qJD(5);
t66 = cos(qJ(5));
t88 = t50 * t57;
t80 = mrSges(6,3) * t88;
t32 = -mrSges(6,2) * t43 + t66 * t80;
t62 = sin(qJ(5));
t81 = qJD(5) * t50;
t34 = (t49 * t62 + t66 * t81) * t57;
t42 = t49 * t60 + qJDD(5);
t76 = t67 * t24 - t25 * t63;
t19 = pkin(4) * t49 + t48 * t90 + t76;
t38 = -t44 * t58 + t55 * t61;
t71 = t19 * t60 + t38 * t57;
t79 = (-mrSges(6,1) * t66 + mrSges(6,2) * t62) * t88 ^ 2;
t15 = m(6) * (-t62 * t20 + t66 * t71) - t34 * mrSges(6,3) + t42 * mrSges(6,1) - t62 * t79 + t43 * t32;
t31 = mrSges(6,1) * t43 - t62 * t80;
t35 = (t49 * t66 - t62 * t81) * t57;
t16 = m(6) * (t66 * t20 + t62 * t71) + t35 * mrSges(6,3) - t42 * mrSges(6,2) + t66 * t79 - t43 * t31;
t91 = t66 * t15 + t62 * t16;
t18 = m(6) * (-t19 * t57 + t38 * t60) + t34 * mrSges(6,2) - t35 * mrSges(6,1) + (t31 * t62 - t32 * t66) * t88;
t74 = m(5) * t38 + t60 * t18 + t57 * t91;
t72 = m(4) * t38 + t74;
t11 = m(3) * t38 + t72;
t12 = m(5) * t83 - t48 * mrSges(5,1) - t49 * mrSges(5,2) - t62 * t15 + t66 * t16;
t9 = m(5) * t76 + t49 * mrSges(5,1) - t48 * mrSges(5,2) - t57 * t18 + t60 * t91;
t7 = m(4) * t75 + t53 * mrSges(4,1) - t52 * mrSges(4,2) + t63 * t12 + t67 * t9;
t8 = m(4) * t82 - t52 * mrSges(4,1) - t53 * mrSges(4,2) + t67 * t12 - t63 * t9;
t5 = m(3) * t73 + qJDD(2) * mrSges(3,1) - t70 * mrSges(3,2) + t64 * t8 + t68 * t7;
t6 = m(3) * t78 - t70 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t64 * t7 + t68 * t8;
t77 = m(2) * t55 + t61 * t11 + t5 * t86 + t6 * t87;
t2 = m(2) * t45 - t5 * t65 + t6 * t69;
t1 = m(2) * t44 - t58 * t11 + (t5 * t69 + t6 * t65) * t61;
t3 = [-m(1) * g(1) - t1 * t56 + t2 * t59, t2, t6, t8, t12, t16; -m(1) * g(2) + t1 * t59 + t2 * t56, t1, t5, t7, t9, t15; -m(1) * g(3) + t77, t77, t11, t72, t74, t18;];
f_new = t3;
