% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2018-11-14 13:55
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4RRRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP1_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:54:21
% EndTime: 2018-11-14 13:54:21
% DurationCPUTime: 0.25s
% Computational Cost: add. (508->60), mult. (1239->83), div. (0->0), fcn. (801->4), ass. (0->43)
t40 = cos(qJ(2));
t36 = pkin(1) * t40 + pkin(2);
t39 = cos(qJ(3));
t37 = sin(qJ(3));
t38 = sin(qJ(2));
t58 = t37 * t38;
t24 = -pkin(1) * t58 + t36 * t39;
t19 = pkin(3) + t24;
t70 = -t19 + t24;
t62 = t39 * pkin(2);
t35 = pkin(3) + t62;
t69 = t35 - t62;
t56 = t38 * t39;
t28 = (-t37 * t40 - t56) * pkin(1);
t26 = t28 * mrSges(5,1);
t27 = t28 * mrSges(4,1);
t68 = (mrSges(3,1) * t38 + mrSges(3,2) * t40) * pkin(1) - t26 - t27;
t63 = t37 * pkin(2);
t67 = mrSges(4,1) + mrSges(5,1);
t53 = mrSges(4,2) + mrSges(5,2);
t65 = m(5) / 0.2e1;
t64 = m(5) * pkin(3);
t25 = pkin(1) * t56 + t36 * t37;
t29 = (t39 * t40 - t58) * pkin(1);
t61 = t25 * t29;
t60 = t29 * mrSges(4,2);
t59 = t29 * mrSges(5,2);
t47 = t53 * t24;
t3 = -t47 + (m(5) * t70 - t67) * t25;
t51 = t3 * qJD(1);
t4 = t53 * t29 - m(4) * (t24 * t28 + t61) - m(5) * (t19 * t28 + t61) + t68;
t50 = t4 * qJD(1);
t49 = -mrSges(4,1) / 0.2e1 - mrSges(5,1) / 0.2e1;
t48 = -mrSges(4,2) / 0.2e1 - mrSges(5,2) / 0.2e1;
t46 = -t67 - t64;
t41 = (-t25 * t69 + t63 * t70) * t65;
t42 = -t27 / 0.2e1 - t26 / 0.2e1 - t28 * t64 / 0.2e1;
t1 = t41 + t42 + t67 * (-t63 / 0.2e1 - t25 / 0.2e1) + t53 * (-t62 / 0.2e1 - t24 / 0.2e1 + t29 / 0.2e1);
t9 = t53 * t62 + (m(5) * t69 + t67) * t63;
t44 = -t1 * qJD(1) + t9 * qJD(2);
t17 = t29 * t63;
t2 = -t60 / 0.2e1 - t59 / 0.2e1 + t49 * t25 + t48 * t24 + (t37 * t49 + t39 * t48) * pkin(2) + t41 - t42;
t5 = [-qJD(2) * t4 + qJD(3) * t3, t2 * qJD(3) - t50 + (-t59 - t60 + m(4) * (t28 * t62 + t17) + 0.2e1 * (t28 * t35 + t17) * t65 - t68) * qJD(2), t2 * qJD(2) + t51 + (t25 * t46 - t47) * qJD(3), 0; qJD(3) * t1 + t50, -t9 * qJD(3) (t37 * t46 - t39 * t53) * qJD(3) * pkin(2) - t44, 0; -qJD(2) * t1 - t51, t44, 0, 0; 0, 0, 0, 0;];
Cq  = t5;
