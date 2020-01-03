% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:36
% EndTime: 2019-12-31 16:33:37
% DurationCPUTime: 0.37s
% Computational Cost: add. (586->71), mult. (1492->103), div. (0->0), fcn. (1229->6), ass. (0->53)
t48 = sin(qJ(4));
t50 = cos(qJ(4));
t94 = -Ifges(5,4) * t48 + (Ifges(5,1) - Ifges(5,2)) * t50;
t96 = t94 * t48;
t47 = t50 ^ 2;
t88 = t48 ^ 2 + t47;
t95 = -t88 * mrSges(5,3) + mrSges(4,2);
t61 = t50 * mrSges(5,1) - t48 * mrSges(5,2);
t91 = mrSges(4,1) + t61;
t79 = cos(qJ(3));
t65 = t79 * pkin(2);
t89 = -t65 / 0.2e1;
t72 = t50 * mrSges(5,2);
t73 = t48 * mrSges(5,1);
t60 = t72 + t73;
t86 = m(5) * pkin(2);
t49 = sin(qJ(3));
t78 = sin(qJ(2));
t80 = cos(qJ(2));
t35 = -t49 * t80 - t79 * t78;
t34 = t49 * t78 - t79 * t80;
t64 = t88 * t34;
t84 = m(5) * (pkin(3) * t35 - pkin(6) * t64);
t83 = t49 * pkin(2);
t17 = (t60 / 0.2e1 + t72 / 0.2e1 + t73 / 0.2e1) * t34;
t5 = m(5) * (-0.1e1 + t88) * t35 * t34;
t70 = t5 * qJD(1);
t81 = t17 * qJD(4) + t70;
t75 = Ifges(5,4) * t50;
t74 = t34 * t49;
t44 = -t65 - pkin(3);
t28 = t44 * t60;
t10 = -Ifges(5,4) * t47 - t28 - t96;
t69 = t10 * qJD(2);
t68 = qJD(2) + qJD(3);
t67 = t84 / 0.2e1;
t63 = Ifges(5,5) * t50 - Ifges(5,6) * t48;
t53 = t88 * t79;
t43 = pkin(6) + t83;
t54 = -t44 * t35 - t43 * t64;
t51 = m(5) * ((-t53 * t35 + t74) * pkin(2) + t54);
t2 = t67 - t51 / 0.2e1;
t52 = -t95 * t65 - t91 * t83;
t9 = (t53 * t43 + t44 * t49) * t86 + t52;
t59 = -t2 * qJD(1) + t9 * qJD(2);
t11 = (pkin(3) * mrSges(5,2) - t75) * t50 + (pkin(3) * mrSges(5,1) - t94) * t48;
t57 = t89 + pkin(3) / 0.2e1;
t7 = -t28 / 0.2e1 + (t57 * mrSges(5,2) - t75) * t50 + (t57 * mrSges(5,1) - t94) * t48;
t58 = -t7 * qJD(2) - t11 * qJD(3);
t56 = t95 * t34 + t91 * t35;
t8 = t75 * t50 + t28 / 0.2e1 + t96 + (-pkin(3) / 0.2e1 + t89) * t60;
t1 = t51 / 0.2e1 + t67 + t56;
t3 = [t68 * t5, (-t80 * mrSges(3,2) - t78 * mrSges(3,1) + m(4) * (t79 * t35 - t74) * pkin(2) + m(5) * t54 + t56) * qJD(2) + t1 * qJD(3) + t81, t1 * qJD(2) + (t56 + t84) * qJD(3) + t81, t61 * qJD(4) * t35 + t68 * t17; -t2 * qJD(3) - t70, t9 * qJD(3) - t10 * qJD(4), ((-pkin(3) * t49 + pkin(6) * t53) * t86 + t52) * qJD(3) + t8 * qJD(4) + t59, -t69 + t8 * qJD(3) + (-t61 * t43 + t63) * qJD(4); t2 * qJD(2) - t70, -t7 * qJD(4) - t59, -t11 * qJD(4), (-t61 * pkin(6) + t63) * qJD(4) + t58; 0, t7 * qJD(3) + t69, -t58, 0;];
Cq = t3;
