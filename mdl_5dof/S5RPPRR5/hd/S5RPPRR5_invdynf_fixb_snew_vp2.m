% Calculate vector of cutting forces with Newton-Euler
% S5RPPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:31
% EndTime: 2019-12-31 17:56:31
% DurationCPUTime: 0.36s
% Computational Cost: add. (3653->85), mult. (5135->109), div. (0->0), fcn. (2168->8), ass. (0->51)
t33 = -qJDD(1) + qJDD(4);
t40 = sin(qJ(5));
t43 = cos(qJ(5));
t34 = -qJD(1) + qJD(4);
t60 = qJD(5) * t34;
t19 = t40 * t33 + t43 * t60;
t20 = t43 * t33 - t40 * t60;
t65 = t34 * t40;
t26 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t65;
t64 = t34 * t43;
t27 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t64;
t32 = t34 ^ 2;
t46 = qJD(1) ^ 2;
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t59 = t42 * g(1) - t45 * g(2);
t24 = qJDD(1) * pkin(1) + t59;
t56 = -t45 * g(1) - t42 * g(2);
t25 = -t46 * pkin(1) + t56;
t38 = sin(pkin(8));
t39 = cos(pkin(8));
t61 = t38 * t24 + t39 * t25;
t57 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t61;
t66 = -pkin(2) - pkin(3);
t14 = t66 * t46 + t57;
t58 = t39 * t24 - t38 * t25;
t49 = -t46 * qJ(3) + qJDD(3) - t58;
t16 = t66 * qJDD(1) + t49;
t41 = sin(qJ(4));
t44 = cos(qJ(4));
t55 = -t41 * t14 + t44 * t16;
t67 = (t40 * t26 - t43 * t27) * t34 + m(6) * (-t33 * pkin(4) - t32 * pkin(7) - t55) - t20 * mrSges(6,1) + t19 * mrSges(6,2);
t63 = mrSges(3,1) + mrSges(4,1);
t62 = t44 * t14 + t41 * t16;
t11 = -t32 * pkin(4) + t33 * pkin(7) + t62;
t18 = (-mrSges(6,1) * t43 + mrSges(6,2) * t40) * t34;
t37 = g(3) - qJDD(2);
t8 = m(6) * (-t40 * t11 + t43 * t37) - t19 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t18 * t65 + qJD(5) * t27;
t9 = m(6) * (t43 * t11 + t40 * t37) + t20 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t18 * t64 - qJD(5) * t26;
t6 = m(5) * t62 - t32 * mrSges(5,1) - t33 * mrSges(5,2) - t40 * t8 + t43 * t9;
t7 = m(5) * t55 + t33 * mrSges(5,1) - t32 * mrSges(5,2) - t67;
t53 = -t41 * t7 + t44 * t6 + m(4) * (-t46 * pkin(2) + t57) + qJDD(1) * mrSges(4,3);
t52 = -m(4) * (-qJDD(1) * pkin(2) + t49) - t41 * t6 - t44 * t7;
t51 = m(5) * t37 + t40 * t9 + t43 * t8;
t48 = -m(4) * t37 - t51;
t47 = -m(3) * t37 + t48;
t4 = m(3) * t58 + (-mrSges(3,2) + mrSges(4,3)) * t46 + t63 * qJDD(1) + t52;
t3 = m(3) * t61 - qJDD(1) * mrSges(3,2) - t63 * t46 + t53;
t2 = m(2) * t56 - t46 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t39 * t3 - t38 * t4;
t1 = m(2) * t59 + qJDD(1) * mrSges(2,1) - t46 * mrSges(2,2) + t38 * t3 + t39 * t4;
t5 = [-m(1) * g(1) - t42 * t1 + t45 * t2, t2, t3, -t46 * mrSges(4,1) + t53, t6, t9; -m(1) * g(2) + t45 * t1 + t42 * t2, t1, t4, t48, t7, t8; (-m(1) - m(2)) * g(3) + t47, -m(2) * g(3) + t47, t47, -qJDD(1) * mrSges(4,1) - t46 * mrSges(4,3) - t52, t51, t67;];
f_new = t5;
