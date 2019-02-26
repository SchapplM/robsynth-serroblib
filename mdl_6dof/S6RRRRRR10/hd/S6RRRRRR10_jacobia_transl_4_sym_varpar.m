% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR10_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobia_transl_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:56
% EndTime: 2019-02-26 22:52:56
% DurationCPUTime: 0.52s
% Computational Cost: add. (356->91), mult. (1016->164), div. (0->0), fcn. (1316->14), ass. (0->63)
t33 = sin(qJ(4));
t37 = cos(qJ(4));
t35 = sin(qJ(2));
t36 = sin(qJ(1));
t39 = cos(qJ(2));
t40 = cos(qJ(1));
t60 = cos(pkin(6));
t56 = t40 * t60;
t21 = t36 * t35 - t39 * t56;
t29 = sin(pkin(7));
t32 = cos(pkin(7));
t30 = sin(pkin(6));
t69 = t30 * t40;
t15 = -t21 * t29 + t32 * t69;
t28 = sin(pkin(8));
t31 = cos(pkin(8));
t22 = t35 * t56 + t36 * t39;
t34 = sin(qJ(3));
t38 = cos(qJ(3));
t59 = t29 * t69;
t5 = (t21 * t32 + t59) * t38 + t22 * t34;
t52 = t15 * t28 + t31 * t5;
t66 = t32 * t34;
t6 = t21 * t66 - t22 * t38 + t34 * t59;
t82 = t52 * t33 + t6 * t37;
t81 = -t6 * t33 + t52 * t37;
t54 = t33 * r_i_i_C(1) + t37 * r_i_i_C(2);
t75 = pkin(12) + r_i_i_C(3);
t58 = t75 * t28;
t78 = -t54 * t31 + t58;
t74 = t29 * pkin(11);
t72 = t28 * t29;
t71 = t29 * t31;
t70 = t30 * t36;
t68 = t31 * t33;
t67 = t31 * t37;
t65 = t32 * t38;
t64 = t34 * t35;
t63 = t34 * t39;
t62 = t35 * t38;
t61 = t38 * t39;
t57 = t36 * t60;
t55 = t60 * t29;
t23 = -t40 * t35 - t39 * t57;
t17 = -t23 * t29 + t32 * t70;
t24 = -t35 * t57 + t40 * t39;
t41 = t23 * t32 + t29 * t70;
t7 = -t24 * t34 + t41 * t38;
t49 = t17 * t28 + t31 * t7;
t13 = t38 * t55 + (t32 * t61 - t64) * t30;
t48 = t13 * t31 + (-t30 * t39 * t29 + t60 * t32) * t28;
t47 = t37 * r_i_i_C(1) - t33 * r_i_i_C(2) + pkin(3);
t9 = t21 * t34 - t22 * t65;
t45 = t22 * t72 + t31 * t9;
t11 = -t23 * t34 - t24 * t65;
t43 = t11 * t31 + t24 * t72;
t14 = t34 * t55 + (t32 * t63 + t62) * t30;
t12 = t23 * t38 - t24 * t66;
t10 = -t21 * t38 - t22 * t66;
t8 = t24 * t38 + t41 * t34;
t2 = t49 * t33 + t8 * t37;
t1 = -t8 * t33 + t49 * t37;
t3 = [t82 * r_i_i_C(1) + t81 * r_i_i_C(2) + t6 * pkin(3) - t22 * pkin(2) - t36 * pkin(1) + pkin(10) * t69 + t15 * pkin(11) + t75 * (t15 * t31 - t5 * t28) (t12 * t37 + t43 * t33) * r_i_i_C(1) + (-t12 * t33 + t43 * t37) * r_i_i_C(2) + t12 * pkin(3) + t23 * pkin(2) + t24 * t74 + t75 * (-t11 * t28 + t24 * t71) (t7 * t37 - t8 * t68) * r_i_i_C(1) + (-t7 * t33 - t8 * t67) * r_i_i_C(2) + t7 * pkin(3) + t8 * t58, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t40 * pkin(1) + t24 * pkin(2) + t8 * pkin(3) + pkin(10) * t70 + t17 * pkin(11) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t75 * (t17 * t31 - t7 * t28) (t10 * t37 + t45 * t33) * r_i_i_C(1) + (-t10 * t33 + t45 * t37) * r_i_i_C(2) + t10 * pkin(3) - t21 * pkin(2) + t22 * t74 + t75 * (t22 * t71 - t9 * t28) (-t37 * t5 + t6 * t68) * r_i_i_C(1) + (t33 * t5 + t6 * t67) * r_i_i_C(2) - t5 * pkin(3) - t6 * t58, -t81 * r_i_i_C(1) + t82 * r_i_i_C(2), 0, 0; 0 (t47 * (-t32 * t64 + t61) - t78 * (-t32 * t62 - t63) + t39 * pkin(2) + (t54 * t28 + t75 * t31 + pkin(11)) * t35 * t29) * t30, t47 * t13 + t78 * t14 (-t14 * t33 + t48 * t37) * r_i_i_C(1) + (-t14 * t37 - t48 * t33) * r_i_i_C(2), 0, 0;];
Ja_transl  = t3;
