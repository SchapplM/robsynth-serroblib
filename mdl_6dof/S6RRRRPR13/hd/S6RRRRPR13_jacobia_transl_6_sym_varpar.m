% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR13_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR13_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:37:34
% EndTime: 2019-02-26 22:37:35
% DurationCPUTime: 0.29s
% Computational Cost: add. (433->72), mult. (1116->122), div. (0->0), fcn. (1459->12), ass. (0->45)
t39 = sin(qJ(2));
t43 = cos(qJ(2));
t55 = cos(pkin(6));
t64 = cos(qJ(1));
t49 = t55 * t64;
t63 = sin(qJ(1));
t28 = t39 * t49 + t43 * t63;
t38 = sin(qJ(3));
t42 = cos(qJ(3));
t35 = sin(pkin(6));
t52 = t35 * t64;
t16 = t28 * t42 - t38 * t52;
t27 = t39 * t63 - t43 * t49;
t37 = sin(qJ(4));
t41 = cos(qJ(4));
t3 = t16 * t37 - t27 * t41;
t4 = t16 * t41 + t27 * t37;
t54 = -r_i_i_C(3) - pkin(11) + pkin(10);
t66 = t42 * pkin(3) + t38 * t54 + pkin(2);
t36 = sin(qJ(6));
t40 = cos(qJ(6));
t65 = pkin(4) + pkin(5);
t46 = r_i_i_C(1) * t40 - r_i_i_C(2) * t36 + t65;
t47 = r_i_i_C(1) * t36 + r_i_i_C(2) * t40 + qJ(5);
t44 = t37 * t47 + t41 * t46 + pkin(3);
t60 = t35 * t39;
t59 = t35 * t43;
t58 = t37 * t42;
t57 = t41 * t42;
t56 = t42 * t43;
t51 = t35 * t63;
t50 = -t28 * t38 - t42 * t52;
t48 = t55 * t63;
t30 = -t39 * t48 + t43 * t64;
t29 = t39 * t64 + t43 * t48;
t26 = t38 * t55 + t42 * t60;
t20 = t30 * t42 + t38 * t51;
t19 = t30 * t38 - t42 * t51;
t14 = t26 * t41 - t37 * t59;
t13 = t26 * t37 + t41 * t59;
t8 = t20 * t41 + t29 * t37;
t7 = t20 * t37 - t29 * t41;
t2 = t36 * t7 + t40 * t8;
t1 = -t36 * t8 + t40 * t7;
t5 = [-t63 * pkin(1) - t28 * pkin(2) - pkin(3) * t16 + pkin(8) * t52 - t27 * pkin(9) - t3 * t47 - t4 * t46 + t54 * t50, t30 * pkin(9) + t47 * (-t29 * t58 - t30 * t41) + t46 * (-t29 * t57 + t30 * t37) - t66 * t29, -t19 * t44 + t20 * t54, -t46 * t7 + t47 * t8, t7, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; pkin(1) * t64 + t30 * pkin(2) + t20 * pkin(3) + pkin(8) * t51 + t29 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(5) + t19 * t54 + t65 * t8, t28 * pkin(9) + t47 * (-t27 * t58 - t28 * t41) + t46 * (-t27 * t57 + t28 * t37) - t66 * t27, t16 * t54 + t44 * t50, -t3 * t46 + t4 * t47, t3 (t3 * t40 - t36 * t4) * r_i_i_C(1) + (-t3 * t36 - t4 * t40) * r_i_i_C(2); 0 (t47 * (t37 * t56 - t39 * t41) + t46 * (t37 * t39 + t41 * t56) + t39 * pkin(9) + t66 * t43) * t35, t54 * t26 + t44 * (-t38 * t60 + t42 * t55) -t13 * t46 + t14 * t47, t13 (t13 * t40 - t14 * t36) * r_i_i_C(1) + (-t13 * t36 - t14 * t40) * r_i_i_C(2);];
Ja_transl  = t5;
