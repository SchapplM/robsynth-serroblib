% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR6
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:05
% EndTime: 2019-02-26 20:07:05
% DurationCPUTime: 0.31s
% Computational Cost: add. (517->88), mult. (1188->159), div. (0->0), fcn. (1549->16), ass. (0->62)
t41 = sin(pkin(7));
t47 = sin(qJ(2));
t49 = cos(qJ(2));
t43 = cos(pkin(12));
t69 = cos(pkin(6));
t63 = t43 * t69;
t67 = sin(pkin(12));
t53 = t47 * t67 - t49 * t63;
t68 = cos(pkin(7));
t42 = sin(pkin(6));
t72 = t42 * t43;
t78 = t41 * t72 + t53 * t68;
t57 = t69 * t67;
t54 = t43 * t47 + t49 * t57;
t64 = t42 * t67;
t77 = -t41 * t64 + t54 * t68;
t76 = r_i_i_C(3) + pkin(11);
t75 = cos(qJ(3));
t39 = pkin(13) + qJ(5);
t37 = sin(t39);
t74 = t37 * t41;
t38 = cos(t39);
t73 = t38 * t41;
t71 = t42 * t47;
t70 = t42 * t49;
t66 = t41 * t71;
t46 = sin(qJ(3));
t62 = t46 * t68;
t61 = t69 * t41;
t59 = (pkin(4) * sin(pkin(13)) + pkin(9)) * t41;
t58 = t68 * t75;
t45 = sin(qJ(6));
t48 = cos(qJ(6));
t56 = r_i_i_C(1) * t48 - r_i_i_C(2) * t45 + pkin(5);
t44 = -pkin(10) - qJ(4);
t55 = r_i_i_C(1) * t45 + r_i_i_C(2) * t48 - t44;
t36 = cos(pkin(13)) * pkin(4) + pkin(3);
t52 = -t37 * t76 - t38 * t56 - t36;
t31 = t43 * t49 - t47 * t57;
t30 = t47 * t63 + t49 * t67;
t29 = -t41 * t70 + t68 * t69;
t28 = (-t47 * t62 + t49 * t75) * t42;
t27 = (t46 * t49 + t47 * t58) * t42;
t24 = t41 * t54 + t64 * t68;
t23 = t41 * t53 - t68 * t72;
t22 = t46 * t61 + (t47 * t75 + t49 * t62) * t42;
t21 = t46 * t71 - t58 * t70 - t61 * t75;
t20 = t28 * t38 + t37 * t66;
t18 = -t31 * t62 - t54 * t75;
t17 = t31 * t58 - t46 * t54;
t16 = -t30 * t62 - t53 * t75;
t15 = t30 * t58 - t46 * t53;
t14 = t31 * t75 - t46 * t77;
t13 = t31 * t46 + t75 * t77;
t12 = t30 * t75 - t46 * t78;
t11 = t30 * t46 + t75 * t78;
t10 = t22 * t38 + t29 * t37;
t8 = t18 * t38 + t31 * t74;
t6 = t16 * t38 + t30 * t74;
t4 = t14 * t38 + t24 * t37;
t2 = t12 * t38 + t23 * t37;
t1 = [0 (t17 * t45 + t48 * t8) * r_i_i_C(1) + (t17 * t48 - t45 * t8) * r_i_i_C(2) + t8 * pkin(5) + t18 * t36 - t17 * t44 - t54 * pkin(2) + t76 * (t18 * t37 - t31 * t73) + t31 * t59, t13 * t52 + t14 * t55, t13, t76 * t4 + t56 * (-t14 * t37 + t24 * t38) (t13 * t48 - t4 * t45) * r_i_i_C(1) + (-t13 * t45 - t4 * t48) * r_i_i_C(2); 0 (t15 * t45 + t48 * t6) * r_i_i_C(1) + (t15 * t48 - t45 * t6) * r_i_i_C(2) + t6 * pkin(5) + t16 * t36 - t15 * t44 - t53 * pkin(2) + t76 * (t16 * t37 - t30 * t73) + t30 * t59, t11 * t52 + t12 * t55, t11, t76 * t2 + t56 * (-t12 * t37 + t23 * t38) (t11 * t48 - t2 * t45) * r_i_i_C(1) + (-t11 * t45 - t2 * t48) * r_i_i_C(2); 1 (t20 * t48 + t27 * t45) * r_i_i_C(1) + (-t20 * t45 + t27 * t48) * r_i_i_C(2) + t20 * pkin(5) + t28 * t36 - t27 * t44 + t76 * (t28 * t37 - t38 * t66) + (t49 * pkin(2) + t47 * t59) * t42, t21 * t52 + t22 * t55, t21, t76 * t10 + t56 * (-t22 * t37 + t29 * t38) (-t10 * t45 + t21 * t48) * r_i_i_C(1) + (-t10 * t48 - t21 * t45) * r_i_i_C(2);];
Ja_transl  = t1;
