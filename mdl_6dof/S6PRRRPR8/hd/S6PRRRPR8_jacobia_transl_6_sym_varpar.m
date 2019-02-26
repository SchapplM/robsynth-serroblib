% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR8
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:36
% EndTime: 2019-02-26 20:14:37
% DurationCPUTime: 0.28s
% Computational Cost: add. (493->73), mult. (1381->138), div. (0->0), fcn. (1805->14), ass. (0->53)
t36 = sin(pkin(7));
t42 = sin(qJ(2));
t45 = cos(qJ(2));
t38 = cos(pkin(12));
t68 = cos(pkin(6));
t61 = t38 * t68;
t66 = sin(pkin(12));
t52 = t66 * t42 - t45 * t61;
t67 = cos(pkin(7));
t37 = sin(pkin(6));
t71 = t37 * t38;
t77 = t36 * t71 + t52 * t67;
t56 = t68 * t66;
t53 = t38 * t42 + t45 * t56;
t62 = t37 * t66;
t76 = -t36 * t62 + t53 * t67;
t75 = pkin(9) * t36;
t74 = cos(qJ(3));
t40 = sin(qJ(4));
t73 = t36 * t40;
t44 = cos(qJ(4));
t72 = t36 * t44;
t70 = t37 * t42;
t69 = t37 * t45;
t65 = r_i_i_C(3) + pkin(11) + pkin(4);
t64 = t36 * t70;
t41 = sin(qJ(3));
t60 = t41 * t67;
t59 = t68 * t36;
t57 = t67 * t74;
t39 = sin(qJ(6));
t43 = cos(qJ(6));
t55 = t39 * r_i_i_C(1) + t43 * r_i_i_C(2) + qJ(5);
t54 = t43 * r_i_i_C(1) - t39 * r_i_i_C(2) + pkin(5) + pkin(10);
t51 = -t36 * t69 + t68 * t67;
t48 = -t55 * t40 - t65 * t44 - pkin(3);
t47 = t52 * t36 - t67 * t71;
t46 = t53 * t36 + t67 * t62;
t31 = t38 * t45 - t42 * t56;
t30 = t42 * t61 + t66 * t45;
t28 = (-t42 * t60 + t74 * t45) * t37;
t24 = t41 * t59 + (t74 * t42 + t45 * t60) * t37;
t23 = t41 * t70 - t57 * t69 - t74 * t59;
t18 = -t31 * t60 - t53 * t74;
t16 = -t30 * t60 - t52 * t74;
t13 = t24 * t40 - t51 * t44;
t12 = t31 * t74 - t76 * t41;
t11 = t31 * t41 + t76 * t74;
t10 = t30 * t74 - t77 * t41;
t9 = t30 * t41 + t77 * t74;
t3 = t12 * t40 - t46 * t44;
t1 = t10 * t40 - t47 * t44;
t2 = [0, t18 * pkin(3) - t53 * pkin(2) + t31 * t75 + t55 * (t18 * t40 - t31 * t72) + t54 * (t31 * t57 - t53 * t41) + t65 * (t18 * t44 + t31 * t73) t11 * t48 + t12 * t54, t55 * (t12 * t44 + t40 * t46) - t65 * t3, t3 (-t11 * t39 + t3 * t43) * r_i_i_C(1) + (-t11 * t43 - t3 * t39) * r_i_i_C(2); 0, t16 * pkin(3) - t52 * pkin(2) + t30 * t75 + t55 * (t16 * t40 - t30 * t72) + t54 * (t30 * t57 - t41 * t52) + t65 * (t16 * t44 + t30 * t73) t10 * t54 + t48 * t9, t55 * (t10 * t44 + t40 * t47) - t65 * t1, t1 (t1 * t43 - t9 * t39) * r_i_i_C(1) + (-t1 * t39 - t9 * t43) * r_i_i_C(2); 1, t28 * pkin(3) + t55 * (t28 * t40 - t44 * t64) + t65 * (t28 * t44 + t40 * t64) + (t45 * pkin(2) + t42 * t75 + t54 * (t41 * t45 + t42 * t57)) * t37, t23 * t48 + t24 * t54, t55 * (t24 * t44 + t40 * t51) - t65 * t13, t13 (t13 * t43 - t23 * t39) * r_i_i_C(1) + (-t13 * t39 - t23 * t43) * r_i_i_C(2);];
Ja_transl  = t2;
