% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:40
% EndTime: 2019-12-31 17:35:41
% DurationCPUTime: 0.34s
% Computational Cost: add. (587->72), mult. (1495->104), div. (0->0), fcn. (1231->6), ass. (0->60)
t84 = cos(qJ(4));
t69 = t84 * pkin(3);
t51 = cos(qJ(5));
t48 = t51 ^ 2;
t49 = sin(qJ(5));
t96 = t49 ^ 2 + t48;
t104 = -t96 * mrSges(6,3) + mrSges(5,2);
t57 = Ifges(6,4) * t49 - (Ifges(6,1) - Ifges(6,2)) * t51;
t102 = -t69 / 0.2e1;
t62 = t51 * mrSges(6,1) - t49 * mrSges(6,2);
t100 = mrSges(5,1) + t62;
t98 = mrSges(6,2) * t102;
t97 = mrSges(6,1) * t102;
t77 = t51 * mrSges(6,2);
t78 = t49 * mrSges(6,1);
t38 = t77 + t78;
t45 = -t69 - pkin(4);
t28 = t45 * t38;
t54 = -Ifges(6,4) * t48 + t57 * t49;
t10 = -t28 + t54;
t50 = sin(qJ(4));
t83 = sin(qJ(3));
t85 = cos(qJ(3));
t34 = t50 * t83 - t84 * t85;
t59 = t77 / 0.2e1 + t78 / 0.2e1;
t16 = (-t38 / 0.2e1 + t59) * t34;
t73 = t16 * qJD(2);
t95 = -t10 * qJD(3) - t73;
t89 = pkin(4) * t38;
t11 = t54 + t89;
t63 = t89 / 0.2e1 - t28 / 0.2e1;
t80 = Ifges(6,4) * t51;
t7 = (-t80 + t98) * t51 + (t57 + t97) * t49 + t63;
t94 = t7 * qJD(3) + t11 * qJD(4) + t73;
t92 = m(6) * pkin(3);
t35 = -t50 * t85 - t84 * t83;
t68 = t96 * t34;
t90 = m(6) * (pkin(4) * t35 - pkin(7) * t68);
t88 = t50 * pkin(3);
t5 = m(6) * (-0.1e1 + t96) * t35 * t34;
t75 = t5 * qJD(2);
t87 = -t16 * qJD(5) - t75;
t17 = (t38 / 0.2e1 + t59) * t34;
t86 = t17 * qJD(5) + t75;
t79 = t34 * t50;
t72 = qJD(3) + qJD(4);
t71 = t90 / 0.2e1;
t67 = Ifges(6,5) * t51 - Ifges(6,6) * t49;
t55 = t96 * t84;
t44 = pkin(7) + t88;
t56 = -t45 * t35 - t44 * t68;
t52 = m(6) * ((-t55 * t35 + t79) * pkin(3) + t56);
t2 = t71 - t52 / 0.2e1;
t53 = -t100 * t88 - t104 * t69;
t9 = (t55 * t44 + t45 * t50) * t92 + t53;
t61 = -t2 * qJD(2) + t9 * qJD(3);
t58 = t100 * t35 + t104 * t34;
t8 = -t63 + (-t57 + t97) * t49 + (t98 + t80) * t51;
t1 = t52 / 0.2e1 + t71 + t58;
t3 = [0, 0, 0, 0, t38 * qJD(5); 0, t72 * t5, (m(6) * t56 + m(5) * (t84 * t35 - t79) * pkin(3) - t85 * mrSges(4,2) - t83 * mrSges(4,1) + t58) * qJD(3) + t1 * qJD(4) + t86, t1 * qJD(3) + (t58 + t90) * qJD(4) + t86, t62 * qJD(5) * t35 + t72 * t17; 0, -t2 * qJD(4) + t87, t9 * qJD(4) - t10 * qJD(5), ((-pkin(4) * t50 + t55 * pkin(7)) * t92 + t53) * qJD(4) + t8 * qJD(5) + t61, t8 * qJD(4) + (-t62 * t44 + t67) * qJD(5) + t95; 0, t2 * qJD(3) + t87, -t7 * qJD(5) - t61, -t11 * qJD(5), (-t62 * pkin(7) + t67) * qJD(5) - t94; 0, t72 * t16, t7 * qJD(4) - t95, t94, 0;];
Cq = t3;
