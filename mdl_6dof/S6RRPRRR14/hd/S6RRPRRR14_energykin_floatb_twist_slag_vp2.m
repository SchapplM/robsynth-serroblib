% Calculate kinetic energy for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S6RRPRRR14_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:09:45
% EndTime: 2018-12-10 18:09:47
% DurationCPUTime: 2.46s
% Computational Cost: add. (16017->171), mult. (25807->260), div. (0->0), fcn. (22854->18), ass. (0->79)
t85 = sin(qJ(1));
t90 = cos(qJ(1));
t66 = t85 * V_base(5) + t90 * V_base(4);
t105 = pkin(10) * t66;
t70 = pkin(9) * V_base(5) + V_base(1);
t71 = -pkin(9) * V_base(4) + V_base(2);
t63 = -t70 * t85 + t90 * t71;
t72 = V_base(6) + qJD(1);
t80 = cos(pkin(6));
t57 = pkin(1) * t72 - t80 * t105 + t63;
t65 = -t85 * V_base(4) + t90 * V_base(5);
t76 = sin(pkin(6));
t61 = -pkin(1) * t65 - t76 * t105 + V_base(3);
t108 = t57 * t80 + t61 * t76;
t64 = t90 * t70 + t85 * t71;
t92 = t65 * t80 + t72 * t76;
t54 = t92 * pkin(10) + t64;
t84 = sin(qJ(2));
t89 = cos(qJ(2));
t45 = t108 * t89 - t54 * t84;
t62 = -t65 * t76 + t72 * t80 + qJD(2);
t79 = cos(pkin(7));
t56 = t66 * t89 + t92 * t84;
t96 = qJ(3) * t56;
t40 = pkin(2) * t62 - t79 * t96 + t45;
t49 = -t57 * t76 + t80 * t61;
t55 = -t66 * t84 + t92 * t89;
t75 = sin(pkin(7));
t44 = -pkin(2) * t55 - t75 * t96 + t49;
t107 = t40 * t79 + t44 * t75;
t73 = sin(pkin(14));
t77 = cos(pkin(14));
t93 = t55 * t79 + t62 * t75;
t48 = t56 * t77 + t93 * t73;
t104 = pkin(11) * t48;
t46 = t108 * t84 + t89 * t54;
t39 = t93 * qJ(3) + t46;
t26 = t107 * t77 - t39 * t73;
t50 = -t55 * t75 + t62 * t79;
t78 = cos(pkin(8));
t22 = pkin(3) * t50 - t78 * t104 + t26;
t31 = -t40 * t75 + t79 * t44 + qJD(3);
t47 = -t56 * t73 + t93 * t77;
t74 = sin(pkin(8));
t25 = -pkin(3) * t47 - t74 * t104 + t31;
t106 = t22 * t78 + t25 * t74;
t27 = t107 * t73 + t77 * t39;
t94 = t47 * t78 + t50 * t74;
t21 = t94 * pkin(11) + t27;
t83 = sin(qJ(4));
t88 = cos(qJ(4));
t11 = t106 * t88 - t83 * t21;
t33 = -t48 * t83 + t94 * t88;
t12 = t106 * t83 + t88 * t21;
t38 = -t47 * t74 + t50 * t78 + qJD(4);
t10 = pkin(12) * t38 + t12;
t15 = -t22 * t74 + t78 * t25;
t34 = t48 * t88 + t94 * t83;
t14 = -pkin(4) * t33 - pkin(12) * t34 + t15;
t82 = sin(qJ(5));
t87 = cos(qJ(5));
t6 = t87 * t10 + t82 * t14;
t5 = -t10 * t82 + t14 * t87;
t29 = -t34 * t82 + t38 * t87;
t9 = -pkin(4) * t38 - t11;
t91 = V_base(3) ^ 2;
t86 = cos(qJ(6));
t81 = sin(qJ(6));
t32 = qJD(5) - t33;
t30 = t34 * t87 + t38 * t82;
t28 = qJD(6) - t29;
t17 = t30 * t86 + t32 * t81;
t16 = -t30 * t81 + t32 * t86;
t7 = -pkin(5) * t29 - pkin(13) * t30 + t9;
t4 = pkin(13) * t32 + t6;
t3 = -pkin(5) * t32 - t5;
t2 = t4 * t86 + t7 * t81;
t1 = -t4 * t81 + t7 * t86;
t8 = (V_base(3) * mrSges(2,2) - t63 * mrSges(2,3) + Ifges(2,5) * t72 + Ifges(2,1) * t66 / 0.2e1) * t66 + (-t31 * mrSges(4,1) + t27 * mrSges(4,3) + Ifges(4,4) * t48 + Ifges(4,6) * t50 + Ifges(4,2) * t47 / 0.2e1) * t47 + (t63 * mrSges(2,1) - t64 * mrSges(2,2) + Ifges(2,3) * t72 / 0.2e1) * t72 + (-t9 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t30 + Ifges(6,6) * t32 + Ifges(6,2) * t29 / 0.2e1) * t29 + (t15 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,5) * t38 + Ifges(5,1) * t34 / 0.2e1) * t34 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t31 * mrSges(4,2) - t26 * mrSges(4,3) + Ifges(4,5) * t50 + Ifges(4,1) * t48 / 0.2e1) * t48 + (t26 * mrSges(4,1) - t27 * mrSges(4,2) + Ifges(4,3) * t50 / 0.2e1) * t50 + (t9 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t32 + Ifges(6,1) * t30 / 0.2e1) * t30 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t28 / 0.2e1) * t28 + (t49 * mrSges(3,2) - t45 * mrSges(3,3) + Ifges(3,5) * t62 + Ifges(3,1) * t56 / 0.2e1) * t56 + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t32 / 0.2e1) * t32 + (-t15 * mrSges(5,1) + t12 * mrSges(5,3) + Ifges(5,4) * t34 + Ifges(5,6) * t38 + Ifges(5,2) * t33 / 0.2e1) * t33 + (-t3 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t17 + Ifges(7,6) * t28 + Ifges(7,2) * t16 / 0.2e1) * t16 + (-V_base(3) * mrSges(2,1) + t64 * mrSges(2,3) + Ifges(2,4) * t66 + Ifges(2,6) * t72 + Ifges(2,2) * t65 / 0.2e1) * t65 + (t3 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t28 + Ifges(7,1) * t17 / 0.2e1) * t17 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t91) / 0.2e1 + m(2) * (t63 ^ 2 + t64 ^ 2 + t91) / 0.2e1 + m(3) * (t45 ^ 2 + t46 ^ 2 + t49 ^ 2) / 0.2e1 + m(4) * (t26 ^ 2 + t27 ^ 2 + t31 ^ 2) / 0.2e1 + (t45 * mrSges(3,1) - t46 * mrSges(3,2) + Ifges(3,3) * t62 / 0.2e1) * t62 + (-t49 * mrSges(3,1) + t46 * mrSges(3,3) + Ifges(3,4) * t56 + Ifges(3,6) * t62 + Ifges(3,2) * t55 / 0.2e1) * t55 + (t11 * mrSges(5,1) - t12 * mrSges(5,2) + Ifges(5,3) * t38 / 0.2e1) * t38 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(6) * (t5 ^ 2 + t6 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t11 ^ 2 + t12 ^ 2 + t15 ^ 2) / 0.2e1;
T  = t8;
