% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPPRRR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPPRRR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobia_transl_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:43
% EndTime: 2019-02-26 19:38:44
% DurationCPUTime: 0.28s
% Computational Cost: add. (756->60), mult. (2166->115), div. (0->0), fcn. (2894->18), ass. (0->66)
t65 = sin(pkin(13));
t66 = sin(pkin(12));
t54 = t66 * t65;
t71 = cos(pkin(13));
t72 = cos(pkin(12));
t63 = t72 * t71;
t75 = cos(pkin(6));
t44 = t75 * t63 - t54;
t74 = cos(pkin(7));
t41 = t44 * t74;
t56 = t66 * t71;
t61 = t72 * t65;
t45 = t75 * t61 + t56;
t68 = sin(pkin(7));
t70 = cos(pkin(14));
t58 = t68 * t70;
t69 = sin(pkin(6));
t50 = t69 * t58;
t64 = sin(pkin(14));
t32 = -t70 * t41 + t45 * t64 + t72 * t50;
t62 = t72 * t69;
t38 = -t44 * t68 - t74 * t62;
t67 = sin(pkin(8));
t73 = cos(pkin(8));
t80 = t32 * t73 - t38 * t67;
t46 = -t75 * t56 - t61;
t42 = t46 * t74;
t47 = -t75 * t54 + t63;
t33 = -t70 * t42 + t47 * t64 - t66 * t50;
t55 = t66 * t69;
t39 = -t46 * t68 + t74 * t55;
t79 = t33 * t73 - t39 * t67;
t60 = t69 * t71;
t51 = t74 * t60;
t59 = t69 * t65;
t37 = -t70 * t51 - t75 * t58 + t64 * t59;
t43 = -t68 * t60 + t75 * t74;
t78 = t37 * t73 - t43 * t67;
t77 = pkin(11) + r_i_i_C(3);
t76 = cos(qJ(4));
t57 = t68 * t64;
t25 = sin(qJ(6));
t28 = cos(qJ(6));
t53 = t28 * r_i_i_C(1) - t25 * r_i_i_C(2) + pkin(5);
t52 = t25 * r_i_i_C(1) + t28 * r_i_i_C(2) + pkin(10);
t49 = t69 * t57;
t26 = sin(qJ(5));
t29 = cos(qJ(5));
t48 = -t77 * t26 - t53 * t29 - pkin(4);
t27 = sin(qJ(4));
t23 = t64 * t51 + t75 * t57 + t70 * t59;
t19 = t37 * t67 + t43 * t73;
t18 = t64 * t42 + t47 * t70 + t66 * t49;
t17 = t64 * t41 + t45 * t70 - t72 * t49;
t14 = t33 * t67 + t39 * t73;
t13 = t32 * t67 + t38 * t73;
t12 = t23 * t76 - t78 * t27;
t11 = t23 * t27 + t78 * t76;
t10 = t18 * t76 - t79 * t27;
t9 = t18 * t27 + t79 * t76;
t8 = t17 * t76 - t80 * t27;
t7 = t17 * t27 + t80 * t76;
t6 = t12 * t29 + t19 * t26;
t4 = t10 * t29 + t14 * t26;
t2 = t13 * t26 + t8 * t29;
t1 = [0, t55, t39, t52 * t10 + t48 * t9, t77 * t4 + t53 * (-t10 * t26 + t14 * t29) (-t4 * t25 + t9 * t28) * r_i_i_C(1) + (-t9 * t25 - t4 * t28) * r_i_i_C(2); 0, -t62, t38, t48 * t7 + t52 * t8, t77 * t2 + t53 * (t13 * t29 - t8 * t26) (-t2 * t25 + t7 * t28) * r_i_i_C(1) + (-t2 * t28 - t7 * t25) * r_i_i_C(2); 1, t75, t43, t48 * t11 + t52 * t12, t77 * t6 + t53 * (-t12 * t26 + t19 * t29) (t11 * t28 - t6 * t25) * r_i_i_C(1) + (-t11 * t25 - t6 * t28) * r_i_i_C(2);];
Ja_transl  = t1;
