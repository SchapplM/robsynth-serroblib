% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR3
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
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:13
% EndTime: 2019-02-26 20:05:14
% DurationCPUTime: 0.35s
% Computational Cost: add. (646->100), mult. (1769->183), div. (0->0), fcn. (2343->16), ass. (0->64)
t48 = sin(pkin(13));
t57 = sin(qJ(3));
t61 = cos(qJ(3));
t70 = cos(pkin(13));
t44 = -t61 * t48 - t57 * t70;
t50 = sin(pkin(7));
t34 = t44 * t50;
t53 = cos(pkin(7));
t36 = t44 * t53;
t51 = sin(pkin(6));
t54 = cos(pkin(6));
t58 = sin(qJ(2));
t62 = cos(qJ(2));
t65 = -t57 * t48 + t61 * t70;
t19 = -t54 * t34 + (-t36 * t62 + t58 * t65) * t51;
t81 = r_i_i_C(3) + pkin(11);
t49 = sin(pkin(12));
t80 = t49 * t51;
t56 = sin(qJ(5));
t79 = t50 * t56;
t60 = cos(qJ(5));
t78 = t50 * t60;
t52 = cos(pkin(12));
t77 = t51 * t52;
t76 = t51 * t53;
t75 = t51 * t58;
t74 = t51 * t62;
t72 = t54 * t58;
t71 = t54 * t62;
t69 = t50 * t75;
t55 = sin(qJ(6));
t59 = cos(qJ(6));
t67 = t59 * r_i_i_C(1) - t55 * r_i_i_C(2) + pkin(5);
t66 = -t55 * r_i_i_C(1) - t59 * r_i_i_C(2) - pkin(10);
t39 = -t49 * t58 + t52 * t71;
t40 = t49 * t62 + t52 * t72;
t64 = -t34 * t77 + t39 * t36 - t40 * t65;
t41 = -t49 * t71 - t52 * t58;
t42 = -t49 * t72 + t52 * t62;
t14 = -t34 * t80 - t41 * t36 + t42 * t65;
t63 = t81 * t56 + t67 * t60 + pkin(4);
t47 = t61 * pkin(3) + pkin(2);
t38 = -t50 * t74 + t54 * t53;
t37 = t53 * t57 * pkin(3) + (-pkin(9) - qJ(4)) * t50;
t35 = t65 * t53;
t33 = t65 * t50;
t29 = -t41 * t50 + t49 * t76;
t28 = -t39 * t50 - t52 * t76;
t27 = (t36 * t58 + t62 * t65) * t51;
t26 = -t35 * t75 + t44 * t74;
t25 = t27 * t60 + t56 * t69;
t23 = t42 * t36 + t41 * t65;
t22 = -t42 * t35 + t41 * t44;
t21 = t40 * t36 + t39 * t65;
t20 = -t40 * t35 + t39 * t44;
t18 = t54 * t33 + (t35 * t62 + t44 * t58) * t51;
t16 = t19 * t60 + t38 * t56;
t13 = t33 * t80 + t41 * t35 + t42 * t44;
t10 = -t33 * t77 + t39 * t35 + t40 * t44;
t8 = t23 * t60 + t42 * t79;
t6 = t21 * t60 + t40 * t79;
t4 = t14 * t60 + t29 * t56;
t2 = t28 * t56 - t60 * t64;
t1 = [0 (-t22 * t55 + t8 * t59) * r_i_i_C(1) + (-t22 * t59 - t8 * t55) * r_i_i_C(2) + t8 * pkin(5) + t23 * pkin(4) - t22 * pkin(10) + t41 * t47 - t42 * t37 + t81 * (t23 * t56 - t42 * t78) -t66 * t14 + (-t42 * t57 + (t41 * t53 + t50 * t80) * t61) * pkin(3) + t63 * t13, t29, t81 * t4 + t67 * (-t14 * t56 + t29 * t60) (-t13 * t59 - t4 * t55) * r_i_i_C(1) + (t13 * t55 - t4 * t59) * r_i_i_C(2); 0 (-t20 * t55 + t6 * t59) * r_i_i_C(1) + (-t20 * t59 - t6 * t55) * r_i_i_C(2) + t6 * pkin(5) + t21 * pkin(4) - t20 * pkin(10) + t39 * t47 - t40 * t37 + t81 * (t21 * t56 - t40 * t78) t66 * t64 + (-t40 * t57 + (t39 * t53 - t50 * t77) * t61) * pkin(3) + t63 * t10, t28, t81 * t2 + t67 * (t28 * t60 + t56 * t64) (-t10 * t59 - t2 * t55) * r_i_i_C(1) + (t10 * t55 - t2 * t59) * r_i_i_C(2); 1 (t25 * t59 - t26 * t55) * r_i_i_C(1) + (-t25 * t55 - t26 * t59) * r_i_i_C(2) + t25 * pkin(5) + t27 * pkin(4) - t26 * pkin(10) + (-t58 * t37 + t62 * t47) * t51 + t81 * (t27 * t56 - t60 * t69) -t66 * t19 + (t54 * t50 * t61 + (t53 * t61 * t62 - t57 * t58) * t51) * pkin(3) + t63 * t18, t38, t81 * t16 + t67 * (-t19 * t56 + t38 * t60) (-t16 * t55 - t18 * t59) * r_i_i_C(1) + (-t16 * t59 + t18 * t55) * r_i_i_C(2);];
Ja_transl  = t1;
