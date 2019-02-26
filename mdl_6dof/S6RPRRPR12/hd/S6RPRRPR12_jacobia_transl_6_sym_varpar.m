% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR12_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR12_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:07:15
% EndTime: 2019-02-26 21:07:16
% DurationCPUTime: 0.33s
% Computational Cost: add. (505->64), mult. (1388->108), div. (0->0), fcn. (1830->14), ass. (0->50)
t38 = cos(qJ(1));
t65 = cos(pkin(12));
t67 = cos(pkin(6));
t55 = t67 * t65;
t63 = sin(pkin(12));
t69 = sin(qJ(1));
t48 = -t38 * t55 + t69 * t63;
t32 = sin(pkin(6));
t64 = sin(pkin(7));
t59 = t32 * t64;
t66 = cos(pkin(7));
t78 = t38 * t59 + t48 * t66;
t53 = t67 * t63;
t25 = t38 * t53 + t69 * t65;
t35 = sin(qJ(3));
t70 = cos(qJ(3));
t14 = -t25 * t70 + t78 * t35;
t34 = sin(qJ(4));
t37 = cos(qJ(4));
t60 = t32 * t66;
t40 = -t38 * t60 + t48 * t64;
t77 = t14 * t34 + t40 * t37;
t76 = t14 * t37 - t40 * t34;
t11 = t25 * t35 + t78 * t70;
t44 = t38 * t63 + t69 * t55;
t72 = t44 * t66 - t69 * t59;
t71 = pkin(5) + pkin(10);
t68 = t38 * t32;
t62 = r_i_i_C(3) + pkin(11) + pkin(4);
t61 = t69 * t32;
t54 = t67 * t64;
t52 = t66 * t65;
t33 = sin(qJ(6));
t36 = cos(qJ(6));
t51 = t33 * r_i_i_C(1) + t36 * r_i_i_C(2) + qJ(5);
t49 = t36 * r_i_i_C(1) - t33 * r_i_i_C(2) + t71;
t43 = -t65 * t59 + t67 * t66;
t42 = -t51 * t34 - t62 * t37 - pkin(3);
t39 = t44 * t64 + t69 * t60;
t26 = t38 * t65 - t69 * t53;
t20 = t35 * t54 + (t35 * t52 + t63 * t70) * t32;
t19 = -t70 * t54 + (t35 * t63 - t52 * t70) * t32;
t16 = t26 * t70 - t72 * t35;
t15 = t26 * t35 + t72 * t70;
t9 = t20 * t34 - t43 * t37;
t8 = t16 * t37 + t39 * t34;
t7 = t16 * t34 - t39 * t37;
t2 = t15 * t36 + t7 * t33;
t1 = -t15 * t33 + t7 * t36;
t3 = [-t69 * pkin(1) - t25 * pkin(2) + t14 * pkin(3) - t40 * pkin(9) + qJ(2) * t68 - t49 * t11 + t51 * t77 + t62 * t76, t61, t42 * t15 + t49 * t16, t51 * t8 - t62 * t7, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t38 * pkin(1) + t26 * pkin(2) + t16 * pkin(3) + t39 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + qJ(2) * t61 + t7 * qJ(5) + t71 * t15 + t62 * t8, -t68, t42 * t11 - t14 * t49, -t51 * t76 + t62 * t77, -t77 (-t11 * t33 - t36 * t77) * r_i_i_C(1) + (-t11 * t36 + t33 * t77) * r_i_i_C(2); 0, t67, t42 * t19 + t49 * t20, -t62 * t9 + t51 * (t20 * t37 + t43 * t34) t9 (-t19 * t33 + t9 * t36) * r_i_i_C(1) + (-t19 * t36 - t9 * t33) * r_i_i_C(2);];
Ja_transl  = t3;
