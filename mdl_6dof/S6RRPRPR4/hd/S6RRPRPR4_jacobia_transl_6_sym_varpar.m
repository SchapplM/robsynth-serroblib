% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:39:32
% EndTime: 2019-02-26 21:39:32
% DurationCPUTime: 0.25s
% Computational Cost: add. (397->70), mult. (792->105), div. (0->0), fcn. (1039->14), ass. (0->47)
t34 = sin(pkin(11));
t40 = sin(qJ(2));
t44 = cos(qJ(2));
t55 = cos(pkin(11));
t24 = -t44 * t34 - t40 * t55;
t41 = sin(qJ(1));
t45 = cos(qJ(1));
t36 = cos(pkin(6));
t48 = -t40 * t34 + t44 * t55;
t47 = t48 * t36;
t10 = t41 * t24 + t45 * t47;
t38 = sin(qJ(6));
t21 = t24 * t36;
t11 = -t45 * t21 + t41 * t48;
t33 = qJ(4) + pkin(12);
t31 = sin(t33);
t32 = cos(t33);
t35 = sin(pkin(6));
t56 = t45 * t35;
t4 = t11 * t32 - t31 * t56;
t42 = cos(qJ(6));
t65 = t10 * t42 + t4 * t38;
t64 = t10 * t38 - t4 * t42;
t43 = cos(qJ(4));
t29 = t43 * pkin(4) + pkin(3);
t51 = t42 * r_i_i_C(1) - t38 * r_i_i_C(2) + pkin(5);
t63 = pkin(10) + r_i_i_C(3);
t46 = t63 * t31 + t51 * t32 + t29;
t62 = t44 * pkin(2);
t59 = t36 * t44;
t57 = t41 * t35;
t39 = sin(qJ(4));
t53 = -t36 * t40 * pkin(2) + (t39 * pkin(4) + pkin(8) + qJ(3)) * t35;
t52 = -t41 * t21 - t45 * t48;
t37 = -qJ(5) - pkin(9);
t50 = t38 * r_i_i_C(1) + t42 * r_i_i_C(2) - t37;
t49 = -t11 * t31 - t32 * t56;
t30 = pkin(1) + t62;
t20 = t24 * t35;
t19 = t48 * t35;
t16 = -t20 * t32 + t36 * t31;
t13 = t45 * t24 - t41 * t47;
t8 = t31 * t57 - t32 * t52;
t7 = -t31 * t52 - t32 * t57;
t2 = -t13 * t38 + t8 * t42;
t1 = -t13 * t42 - t8 * t38;
t3 = [-t4 * pkin(5) + t64 * r_i_i_C(1) + t65 * r_i_i_C(2) - t10 * t37 - t11 * t29 - t41 * t30 + t53 * t45 + t63 * t49 (-t45 * t40 - t41 * t59) * pkin(2) - t50 * t52 + t46 * t13, t57, t63 * t8 + (t39 * t52 + t43 * t57) * pkin(4) - t51 * t7, -t13, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t8 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t13 * t37 - t29 * t52 + t45 * t30 + t53 * t41 + t63 * t7 (-t41 * t40 + t45 * t59) * pkin(2) + t50 * t11 + t46 * t10, -t56, t63 * t4 + (-t11 * t39 - t43 * t56) * pkin(4) + t51 * t49, -t10, -t65 * r_i_i_C(1) + t64 * r_i_i_C(2); 0, t46 * t19 - t50 * t20 + t35 * t62, t36, t63 * t16 + (t20 * t39 + t36 * t43) * pkin(4) + t51 * (t20 * t31 + t36 * t32) -t19 (-t16 * t38 - t19 * t42) * r_i_i_C(1) + (-t16 * t42 + t19 * t38) * r_i_i_C(2);];
Ja_transl  = t3;
