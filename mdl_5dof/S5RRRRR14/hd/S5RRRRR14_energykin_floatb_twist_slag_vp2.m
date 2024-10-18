% Calculate kinetic energy for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR14_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR14_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_energykin_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR14_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR14_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:17
% EndTime: 2024-09-27 18:42:17
% DurationCPUTime: 0.44s
% Computational Cost: add. (3591->137), mult. (5482->205), div. (0->0), fcn. (4556->12), ass. (0->57)
t53 = V_base(5) * pkin(6) + V_base(1);
t54 = -V_base(4) * pkin(6) + V_base(2);
t63 = sin(qJ(1));
t68 = cos(qJ(1));
t44 = -t53 * t63 + t68 * t54;
t48 = t63 * V_base(5) + t68 * V_base(4);
t56 = V_base(6) + qJD(1);
t39 = pkin(1) * t56 - pkin(7) * t48 + t44;
t45 = t68 * t53 + t63 * t54;
t47 = -t63 * V_base(4) + t68 * V_base(5);
t41 = pkin(7) * t47 + t45;
t62 = sin(qJ(2));
t67 = cos(qJ(2));
t32 = t67 * t39 - t41 * t62;
t55 = qJD(2) + t56;
t58 = cos(pkin(5));
t43 = t47 * t62 + t48 * t67;
t73 = pkin(8) * t43;
t25 = pkin(2) * t55 - t58 * t73 + t32;
t42 = t47 * t67 - t48 * t62;
t46 = -pkin(1) * t47 + V_base(3);
t57 = sin(pkin(5));
t29 = -pkin(2) * t42 - t57 * t73 + t46;
t74 = t25 * t58 + t29 * t57;
t33 = t62 * t39 + t67 * t41;
t70 = t42 * t58 + t55 * t57;
t24 = t70 * pkin(8) + t33;
t61 = sin(qJ(3));
t66 = cos(qJ(3));
t16 = t66 * t24 + t74 * t61;
t30 = -t43 * t61 + t70 * t66;
t12 = pkin(9) * t30 + t16;
t60 = sin(qJ(4));
t65 = cos(qJ(4));
t15 = -t24 * t61 + t74 * t66;
t31 = t43 * t66 + t70 * t61;
t36 = -t42 * t57 + t58 * t55 + qJD(3);
t9 = pkin(3) * t36 - pkin(9) * t31 + t15;
t6 = t65 * t12 + t60 * t9;
t5 = -t12 * t60 + t65 * t9;
t18 = -t25 * t57 + t58 * t29;
t17 = -pkin(3) * t30 + t18;
t35 = qJD(4) + t36;
t69 = V_base(3) ^ 2;
t64 = cos(qJ(5));
t59 = sin(qJ(5));
t34 = qJD(5) + t35;
t20 = t30 * t60 + t31 * t65;
t19 = t30 * t65 - t31 * t60;
t14 = t19 * t59 + t20 * t64;
t13 = t19 * t64 - t20 * t59;
t10 = -pkin(4) * t19 + t17;
t4 = pkin(10) * t19 + t6;
t3 = pkin(4) * t35 - pkin(10) * t20 + t5;
t2 = t3 * t59 + t4 * t64;
t1 = t3 * t64 - t4 * t59;
t7 = m(2) * (t44 ^ 2 + t45 ^ 2 + t69) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t69) / 0.2e1 + m(3) * (t32 ^ 2 + t33 ^ 2 + t46 ^ 2) / 0.2e1 + m(5) * (t17 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t18 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t44 * mrSges(2,1) - t45 * mrSges(2,2) + Ifges(2,3) * t56 / 0.2e1) * t56 + (t32 * mrSges(3,1) - t33 * mrSges(3,2) + Ifges(3,3) * t55 / 0.2e1) * t55 + (t15 * mrSges(4,1) - t16 * mrSges(4,2) + Ifges(4,3) * t36 / 0.2e1) * t36 + (t5 * mrSges(5,1) - t6 * mrSges(5,2) + Ifges(5,3) * t35 / 0.2e1) * t35 + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t34 / 0.2e1) * t34 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t44 * mrSges(2,3) + Ifges(2,5) * t56 + Ifges(2,1) * t48 / 0.2e1) * t48 + (t46 * mrSges(3,2) - t32 * mrSges(3,3) + Ifges(3,5) * t55 + Ifges(3,1) * t43 / 0.2e1) * t43 + (t18 * mrSges(4,2) - t15 * mrSges(4,3) + Ifges(4,5) * t36 + Ifges(4,1) * t31 / 0.2e1) * t31 + (t17 * mrSges(5,2) - t5 * mrSges(5,3) + Ifges(5,5) * t35 + Ifges(5,1) * t20 / 0.2e1) * t20 + (t10 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t34 + Ifges(6,1) * t14 / 0.2e1) * t14 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t45 * mrSges(2,3) + Ifges(2,4) * t48 + Ifges(2,6) * t56 + Ifges(2,2) * t47 / 0.2e1) * t47 + (-t46 * mrSges(3,1) + t33 * mrSges(3,3) + Ifges(3,4) * t43 + Ifges(3,6) * t55 + Ifges(3,2) * t42 / 0.2e1) * t42 + (-t18 * mrSges(4,1) + t16 * mrSges(4,3) + Ifges(4,4) * t31 + Ifges(4,6) * t36 + Ifges(4,2) * t30 / 0.2e1) * t30 + (-t17 * mrSges(5,1) + t6 * mrSges(5,3) + Ifges(5,4) * t20 + Ifges(5,6) * t35 + Ifges(5,2) * t19 / 0.2e1) * t19 + (-t10 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t14 + Ifges(6,6) * t34 + Ifges(6,2) * t13 / 0.2e1) * t13;
T = t7;
