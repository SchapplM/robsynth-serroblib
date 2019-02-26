% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:50:41
% EndTime: 2019-02-26 22:50:41
% DurationCPUTime: 0.22s
% Computational Cost: add. (423->62), mult. (623->97), div. (0->0), fcn. (782->14), ass. (0->43)
t37 = sin(qJ(3));
t40 = cos(qJ(3));
t35 = qJ(4) + qJ(5);
t31 = cos(t35);
t23 = pkin(5) * t31 + cos(qJ(4)) * pkin(4);
t21 = pkin(3) + t23;
t32 = qJ(6) + t35;
t28 = sin(t32);
t29 = cos(t32);
t45 = r_i_i_C(1) * t29 - r_i_i_C(2) * t28 + t21;
t54 = r_i_i_C(3) + pkin(12) + pkin(11) + pkin(10);
t59 = t54 * t37 + t45 * t40 + pkin(2);
t38 = sin(qJ(2));
t39 = sin(qJ(1));
t41 = cos(qJ(2));
t49 = cos(pkin(6));
t53 = cos(qJ(1));
t46 = t49 * t53;
t18 = t38 * t46 + t39 * t41;
t36 = sin(pkin(6));
t48 = t36 * t53;
t10 = t18 * t40 - t37 * t48;
t17 = t39 * t38 - t41 * t46;
t58 = (-t10 * t28 + t17 * t29) * r_i_i_C(1) + (-t10 * t29 - t17 * t28) * r_i_i_C(2);
t47 = t39 * t49;
t20 = -t38 * t47 + t53 * t41;
t52 = t36 * t39;
t14 = t20 * t40 + t37 * t52;
t19 = t53 * t38 + t41 * t47;
t5 = -t14 * t28 + t19 * t29;
t6 = t14 * t29 + t19 * t28;
t57 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t51 = t36 * t40;
t16 = t49 * t37 + t38 * t51;
t50 = t36 * t41;
t56 = (-t16 * t28 - t29 * t50) * r_i_i_C(1) + (-t16 * t29 + t28 * t50) * r_i_i_C(2);
t30 = sin(t35);
t22 = pkin(5) * t30 + sin(qJ(4)) * pkin(4);
t55 = pkin(9) + t22;
t44 = t28 * r_i_i_C(1) + t29 * r_i_i_C(2) + t55;
t43 = -t18 * t37 - t40 * t48;
t13 = t20 * t37 - t39 * t51;
t1 = [-t39 * pkin(1) - t18 * pkin(2) + pkin(8) * t48 - t45 * t10 - t44 * t17 + t54 * t43, -t19 * t59 + t44 * t20, -t45 * t13 + t54 * t14, -t14 * t22 + t19 * t23 + t57 (-t14 * t30 + t19 * t31) * pkin(5) + t57, t57; t53 * pkin(1) + t20 * pkin(2) + pkin(8) * t52 + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t54 * t13 + t14 * t21 + t55 * t19, -t17 * t59 + t44 * t18, t54 * t10 + t45 * t43, -t10 * t22 + t17 * t23 + t58 (-t10 * t30 + t17 * t31) * pkin(5) + t58, t58; 0 (t44 * t38 + t59 * t41) * t36, t54 * t16 + t45 * (-t36 * t38 * t37 + t49 * t40) -t16 * t22 - t23 * t50 + t56 (-t16 * t30 - t31 * t50) * pkin(5) + t56, t56;];
Ja_transl  = t1;
