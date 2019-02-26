% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% Ja_transl [3x7]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S7RRRRRRR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobia_transl_6_sym_varpar: qJ has to be [7x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S7RRRRRRR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobia_transl_6_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:20
% EndTime: 2019-02-26 22:54:20
% DurationCPUTime: 0.49s
% Computational Cost: add. (255->98), mult. (711->189), div. (0->0), fcn. (919->12), ass. (0->65)
t48 = cos(qJ(3));
t42 = sin(qJ(3));
t50 = cos(qJ(1));
t55 = t50 * t42;
t44 = sin(qJ(1));
t49 = cos(qJ(2));
t60 = t44 * t49;
t30 = t48 * t60 + t55;
t41 = sin(qJ(4));
t43 = sin(qJ(2));
t47 = cos(qJ(4));
t63 = t43 * t47;
t14 = t30 * t41 - t44 * t63;
t39 = sin(qJ(6));
t64 = t43 * t44;
t15 = t30 * t47 + t41 * t64;
t54 = t50 * t48;
t29 = t42 * t60 - t54;
t40 = sin(qJ(5));
t46 = cos(qJ(5));
t4 = t15 * t46 - t29 * t40;
t45 = cos(qJ(6));
t79 = -t14 * t45 + t4 * t39;
t78 = -t14 * t39 - t4 * t45;
t75 = -t15 * t40 - t29 * t46;
t74 = t40 * r_i_i_C(3);
t73 = t41 * pkin(3);
t70 = t39 * t41;
t69 = t39 * t46;
t68 = t40 * t47;
t67 = t41 * t45;
t66 = t42 * t43;
t65 = t42 * t49;
t62 = t43 * t48;
t61 = t43 * t50;
t59 = t45 * t46;
t58 = t46 * t47;
t57 = t49 * t41;
t56 = t49 * t47;
t53 = t40 * t66;
t52 = t46 * t66;
t51 = t45 * r_i_i_C(1) - t39 * r_i_i_C(2);
t28 = t47 * t62 - t57;
t27 = t41 * t62 + t56;
t34 = -t44 * t42 + t49 * t54;
t33 = -t44 * t48 - t49 * t55;
t32 = t43 * t41 + t48 * t56;
t31 = t48 * t57 - t63;
t25 = t28 * t50;
t24 = t27 * t50;
t23 = t28 * t44;
t22 = t27 * t44;
t20 = t34 * t47 + t41 * t61;
t19 = t34 * t41 - t47 * t61;
t18 = t32 * t46 - t40 * t65;
t13 = t28 * t46 - t53;
t11 = -t25 * t46 + t50 * t53;
t10 = -t23 * t46 + t44 * t53;
t9 = t33 * t58 - t34 * t40;
t8 = -t29 * t58 - t30 * t40;
t7 = t20 * t46 + t33 * t40;
t6 = t20 * t40 - t33 * t46;
t2 = t19 * t39 + t7 * t45;
t1 = t19 * t45 - t7 * t39;
t3 = [pkin(2) * t64 - t14 * pkin(3) + t78 * r_i_i_C(1) + t79 * r_i_i_C(2) + t75 * r_i_i_C(3) (t11 * t45 - t24 * t39) * r_i_i_C(1) + (-t11 * t39 - t24 * t45) * r_i_i_C(2) + (-t25 * t40 - t50 * t52) * r_i_i_C(3) - t24 * pkin(3) - t50 * t49 * pkin(2) (t33 * t70 + t9 * t45) * r_i_i_C(1) + (t33 * t67 - t9 * t39) * r_i_i_C(2) + (t33 * t68 + t34 * t46) * r_i_i_C(3) + t33 * t73 (-t19 * t59 + t20 * t39) * r_i_i_C(1) + (t19 * t69 + t20 * t45) * r_i_i_C(2) - t19 * t74 + t20 * pkin(3), t7 * r_i_i_C(3) - t51 * t6, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; -pkin(2) * t61 + t19 * pkin(3) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t6 * r_i_i_C(3) (t10 * t45 - t22 * t39) * r_i_i_C(1) + (-t10 * t39 - t22 * t45) * r_i_i_C(2) + (-t23 * t40 - t44 * t52) * r_i_i_C(3) - t22 * pkin(3) - pkin(2) * t60 (-t29 * t70 + t8 * t45) * r_i_i_C(1) + (-t29 * t67 - t8 * t39) * r_i_i_C(2) + (-t29 * t68 + t30 * t46) * r_i_i_C(3) - t29 * t73 (-t14 * t59 + t15 * t39) * r_i_i_C(1) + (t14 * t69 + t15 * t45) * r_i_i_C(2) - t14 * t74 + t15 * pkin(3), t4 * r_i_i_C(3) + t51 * t75, -t79 * r_i_i_C(1) + t78 * r_i_i_C(2), 0; 0 (t18 * t45 + t31 * t39) * r_i_i_C(1) + (-t18 * t39 + t31 * t45) * r_i_i_C(2) + (t32 * t40 + t46 * t65) * r_i_i_C(3) + t31 * pkin(3) - t43 * pkin(2) (t51 * (-t40 * t48 - t42 * t58) + (-t42 * t68 + t46 * t48) * r_i_i_C(3) + (-t39 * r_i_i_C(1) - t45 * r_i_i_C(2) - pkin(3)) * t42 * t41) * t43 (-t27 * t59 + t28 * t39) * r_i_i_C(1) + (t27 * t69 + t28 * t45) * r_i_i_C(2) - t27 * t74 + t28 * pkin(3), t13 * r_i_i_C(3) + t51 * (-t28 * t40 - t52) (-t13 * t39 + t27 * t45) * r_i_i_C(1) + (-t13 * t45 - t27 * t39) * r_i_i_C(2), 0;];
Ja_transl  = t3;
